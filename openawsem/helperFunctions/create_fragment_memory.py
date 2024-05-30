import argparse
from pathlib import Path
import sys
import os
from .IndexPdb import *
from .Pdb2GroLib import *
import subprocess
import argparse
from pathlib import Path
import logging
import requests
import gzip
from concurrent.futures import ThreadPoolExecutor, as_completed
import tempfile


def update_failed_pdb_list(failed_pdb, failed_pdb_list_file):
    try:
        with open(failed_pdb_list_file, 'r') as file:
            existing_pdb_ids = sorted(line.strip() for line in file)
    except FileNotFoundError:
        existing_pdb_ids = []

    # Collect failed pdbIDs that are marked as True
    new_pdb_ids = sorted(pdbID for pdbID, failed in failed_pdb.items() if failed)

    # Merge two sorted lists while removing duplicates
    merged_pdb_ids = []
    i, j = 0, 0
    while i < len(existing_pdb_ids) and j < len(new_pdb_ids):
        if existing_pdb_ids[i] == new_pdb_ids[j]:
            merged_pdb_ids.append(existing_pdb_ids[i])
            i += 1
            j += 1
        elif existing_pdb_ids[i] < new_pdb_ids[j]:
            merged_pdb_ids.append(existing_pdb_ids[i])
            i += 1
        else:
            merged_pdb_ids.append(new_pdb_ids[j])
            j += 1

    # Append remaining items if any list was not exhausted
    if i < len(existing_pdb_ids):
        merged_pdb_ids.extend(existing_pdb_ids[i:])
    if j < len(new_pdb_ids):
        merged_pdb_ids.extend(new_pdb_ids[j:])

    # Write the merged list back to the file
    with open(failed_pdb_list_file, 'w') as file:
        for pdbID in merged_pdb_ids:
            file.write(f"{pdbID}\n")

def process_window(i, record, fragment_length, evalue_threshold, database, residue_base):
    rangeStart = i - 1
    rangeEnd = i + fragment_length - 1
    logging.debug(f"window position::: {i}")
    subrange = str(record[rangeStart:rangeEnd].seq)
    logging.debug(f"fragment subrange::: {subrange}")

    # Use a fixed file name as specified

    with tempfile.NamedTemporaryFile(mode='w', delete=True) as temp_fragment_file:
        temp_fragment_file.write(">query\n")
        temp_fragment_file.write(subrange)
        temp_fragment_file.flush()
        
        # Constructing the PSI-BLAST command
        exeline = f"psiblast -num_iterations 5 -comp_based_stats 0 -word_size 2 -evalue {evalue_threshold} " \
                f"-outfmt '6 sseqid qstart qend sstart send qseq sseq length gaps bitscore evalue' -matrix BLOSUM62 -threshold 9 -window_size 0 " \
                f"-db {database} -query {temp_fragment_file.name}"
        logging.debug(f"executing::: {exeline}")

        psiblastOut = os.popen(exeline).read().splitlines()
    
    if psiblastOut and psiblastOut[-1] == 'Search has CONVERGED!':
        psiblastOut = psiblastOut[:-2]  # exclude last two lines for the text

    # Filter last iteration. 
    # If there are multiple iterations, it can be catched if the e_score decreases and the pdbID changes
    N_blast = len(psiblastOut)
    N_start_new = 1
    line_count = 0
    e_score_old = 0
    pdbID_old = 'BBBBB'   # set initial variables for processing PSI-BLAST output
    for line in psiblastOut:  # For PSI-BLAST with multiple Iterations, find N_start_new for the starting position of the output alignments of the final round
        line_count += 1
        if line_count >= N_blast:
            break
        that = line.split()
        pdbID = that[0]
        e_score = float(that[10])
        if e_score < e_score_old:
            N_start_new = line_count
        if pdbID != pdbID_old:
            e_score_old = e_score
            pdbID_old = pdbID
    logging.debug(f"Number of searched PDBs: {N_blast - N_start_new + 1}")
    
    # Convert psiblastOut to a list, sorted by evalue
    psilist = [None] * (N_blast - N_start_new + 1)
    line_count = 0
    kk = 0
    for line in psiblastOut:
        line_count += 1
        if line_count < N_start_new:
            continue
        if line_count > N_blast:
            break
        that = line.split()
        list_tmp = list()
        for ii in range(0, 11):  # PSI-BLAST output has 11 columns
            if not ii == 10:
                list_tmp.append(that[ii])  # column 10 is evalue
            else:
                list_tmp.append(float(that[ii]))
        psilist[kk] = list_tmp
        kk += 1
        logging.debug(list_tmp[:20])
    psilist.sort(key=lambda x: x[10])


    result = []
    # write output alignments to match file
    for jj in range(0, N_blast - N_start_new + 1):
        this = psilist[jj]
        this[10] = str(this[10])
        this.append(str(i))
        queryStart = int(this[1]) + rangeStart + residue_base
        queryEnd = int(this[2]) + rangeStart + residue_base
        this[1] = str(queryStart)
        this[2] = str(queryEnd)
        out = ' '.join(this)
        gaps = this[8]
        if(gaps == '0'):
            result+=[out]  # skip gapped alignments
    
    return result


def download_pdb(pdbID, pdb_dir, max_retries=5):
    pdbID_lower = pdbID.lower()
    filename = f"{pdbID_lower.upper()}.pdb"
    filepath = Path(pdb_dir) / filename

    if filepath.exists() and filepath.stat().st_size > 0:
        logging.info(f"File {filename} already exists, skipping download.")
        return True

    download_url = f"https://files.wwpdb.org/pub/pdb/data/structures/divided/pdb/{pdbID_lower[1:3]}/pdb{pdbID_lower}.ent.gz"
    temp_gz_path = Path(pdb_dir) / f"{pdbID_lower}.ent.gz"

    retry_count = 0
    while retry_count < max_retries:
        try:
            # Download the file with timeout
            response = requests.get(download_url, stream=True, timeout=10)
            response.raise_for_status()  # Check for HTTP errors

            # Write to temporary file
            with open(temp_gz_path, 'wb') as temp_file:
                for chunk in response.iter_content(chunk_size=8192):
                    temp_file.write(chunk)

            # Unzip the file
            with gzip.open(temp_gz_path, 'rb') as gz_file:
                with open(filepath, 'wb') as output_file:
                    output_file.write(gz_file.read())

            logging.info(f"Successfully downloaded and saved {filename} to {filepath}")
            return True

        except requests.exceptions.HTTPError as e:
            logging.warning(f"HTTP error on retry {retry_count + 1} for {pdbID}: {e}")
        except requests.exceptions.ConnectionError as e:
            logging.warning(f"Connection error on retry {retry_count + 1} for {pdbID}: {e}")
        except requests.exceptions.Timeout as e:
            logging.warning(f"Timeout occurred on retry {retry_count + 1} for {pdbID}: {e}")
        except requests.exceptions.RequestException as e:
            logging.warning(f"Request error on retry {retry_count + 1} for {pdbID}: {e}")
        except IOError as e:
            logging.warning(f"File I/O error on retry {retry_count + 1} for {pdbID}: {e}")
        except Exception as e:
            logging.error(f"An unexpected error occurred on retry {retry_count + 1} for {pdbID}: {e}")
        finally:
            retry_count += 1
            if temp_gz_path.exists():
                temp_gz_path.unlink()

    logging.error(f"Failed to download {filename} after {max_retries} retries.")
    return None

def download_pdbs(pdb_ids, pdb_dir):
    failed_pdb = {}

    with ThreadPoolExecutor(max_workers=10) as executor:
        future_to_pdb = {executor.submit(download_pdb, pdb_id, pdb_dir): pdb_id for pdb_id in pdb_ids}
        
        for future in as_completed(future_to_pdb):
            pdb_id = future_to_pdb[future]
            try:
                result = future.result()
                if result:
                    logging.info(f"Download succeeded for {pdb_id}")
                    failed_pdb[pdb_id] = False
                else:
                    logging.info(f"Download failed for {pdb_id}")
                    failed_pdb[pdb_id] = True
            except Exception as e:
                logging.error(f"Error processing {pdb_id}: {str(e)}")
                failed_pdb[pdb_id] = True

    return failed_pdb

def download_pdb_seqres(pdb_seqres):
    logging.debug(f"Checking if {pdb_seqres} exists")
    pdb_seqres=Path(pdb_seqres)
    if not pdb_seqres.exists():
        import urllib
        logging.warning("pdb_seqres.txt was not found. Attempting download from ftp://ftp.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt")
        url = "ftp://ftp.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt"
        try:
            urllib.request.urlretrieve(url, pdb_seqres)
            logging.info(f"Download complete. Saved to {pdb_seqres}")
            return pdb_seqres
        except urllib.error.URLError as e:
            logging.error(f"Error downloading file: {e.reason}")
        except Exception as e:
            logging.error(f"An error occurred: {e}")
        raise FileNotFoundError(f"Failed to download {pdb_seqres}. Make sure it exists in {pdb_seqres} or provide a valid path.")

def create_index_files(iter, line, N_mem, brain_damage,count, failed_pdb,homo, homo_count, weight, frag_lib_dir, pdb_dir, index_dir, pdb_seqres):
    Missing_count=0
    Missing_pdb={}
    out=''
    out1=''
    if not(iter == 0):
        entries = line.split()
        windows_index_str = entries[11]
        if count[windows_index_str] >= N_mem:
            return out, out1, Missing_pdb, Missing_count
        # pdbfull = str(entries[0])
        entry = entries[0]
        if entry[:3] == "pdb":
            # for example 'pdb|4V12|A'
            pdbfull = str(entry[4:8]) + str(entry[9:])
        else:
            pdbfull = str(entry)
        pdbID = pdbfull[0:4].lower()
        chainID = pdbfull[4:5]
        groFile = frag_lib_dir + pdbID + chainID + ".gro"
        groName = pdbID + chainID + ".gro"
        pdbFile = pdb_dir + pdbID.upper() + ".pdb"
        indexFile = index_dir + pdbID + chainID + ".index"

        if failed_pdb[pdbID]:
            return out, out1, Missing_pdb, Missing_count  # failed-downloaded ones are still in matchlines, need to be ignored
        if brain_damage == 2:
            if homo[pdbID]:
                homo_count[pdbID] += 1
                pass
            else:
                print(pdbID, "is not a homolog, discard")
                return out, out1, Missing_pdb, Missing_count

        if homo[pdbID]:
            if brain_damage == 0:
                print(pdbID, " Using  a homolog.")
                homo_count[pdbID] += 1
            if brain_damage == 1:
                print(pdbID, " is a homolog, discard")
                return out, out1, Missing_pdb, Missing_count
        residue_list = entries[6]  # sseq

        res_Start = int(entries[3])
        res_End = int(entries[4])
        print(pdbFile, "start: ", res_Start, "end: ", res_End)
        # Do I have the index file?  If No, write it

        if not os.path.isfile(indexFile):
            # generate fasta file
            fastaFile = pdbID + '_' + chainID
            with tempfile.NamedTemporaryFile(mode='w', delete=True) as temp_file:
                exeline = "grep -A1 " + fastaFile + " " + pdb_seqres + " > " + temp_file.name
                print("generating fastaFile: ", fastaFile)
                # p = os.popen(exeline)
                subprocess.Popen(exeline, shell=True).wait()
                # p_status = p.wait()
                if os.path.getsize(temp_file.name) > 0:
                    writeIndexFile(temp_file.name, pdbFile,
                                indexFile, chainID)
                    print("Writing indexFile: ", indexFile)
        else:
            print(indexFile, "exist, no need to create.")

        if not os.path.isfile(indexFile):
            print("Can't create index file, ignore and go on!")
            return out, out1, Missing_pdb, Missing_count

        # Read index file
        index = open(indexFile, 'r')
        # create new_index for frag_seq starting position
        line_count = 0
        flag = ' '
        index_shift = 0
        # read and get the flag
        indexlines = list()
        for index_line in index.readlines():
            indexlines.append(index_line)
            line_count += 1
            tmp_line = index_line.split()
            if line_count == 1:  # first line is the flag
                flag = tmp_line[0]  # index_line
            if flag == "SHIFT" and line_count == 2:
                index_shift = int(tmp_line[0])
                print("shift: ", tmp_line[0])

        r_list = ''  # list()
        if flag == "SKIP":
            Missing_pdb[pdbID] = 1
            Missing_count += 1
            print("***********", flag)
            print("SKIP pdb:", pdbID + chainID)
            return out, out1, Missing_pdb, Missing_count
        elif flag == "FULLMATCH":
            new_index = int(entries[3])
            r_list = residue_list
            print("***********", flag)
        elif flag == "SHIFT":
            new_index = int(entries[3]) + index_shift
            r_list = residue_list
            print("***********", flag)
        elif flag == "INDEXED":
            print("***********", flag)
            # check if there is gaps for the fragment of sequence
            count_flag = 0
            line_count1 = 0
            for index_line in indexlines:
                line_count1 += 1
                if not line_count1 == 1:
                    index_entries = index_line.split()
                    seq_id = int(index_entries[0])
                    res_id = int(index_entries[1])
                    if seq_id < res_Start:
                        continue
                    if seq_id > res_End:
                        break
                    if res_id == -1:
                        print("Missing residues in PDB: ", pdbID + chainID)
                        break
                    if count_flag == 0:
                        new_index = res_id
                        count_flag += 1
                    res_nm = index_entries[2]
                    r_list += res_nm
        else:
            print("Skip wrongly written index file ", indexFile)
            return out, out1, Missing_pdb, Missing_count

        if r_list != residue_list:
            print("Missing residues: ", pdbID + chainID, residue_list, " incomplete: ", r_list)
            Missing_pdb[pdbID] = 1
            Missing_count += 1
            return out, out1, Missing_pdb, Missing_count

        if os.path.isfile(pdbFile):
            if not os.path.isfile(groFile):
                print("converting...... " + pdbFile + " --> " + groFile + " chainID: " + chainID)
                Pdb2Gro(pdbFile, groFile, chainID)
            else:
                print("Exist " + groFile)
            count[windows_index_str] += 1

            length = res_End - res_Start + 1
            out = groFile + ' ' + entries[1] + ' '  # queue start
            out += str(new_index) + ' ' + str(length) + ' ' + \
                str(weight) + "\n"  # frag_seq start
            out1 = windows_index_str + ' ' + str(count[windows_index_str])
            out1 += ' ' + entries[9] + ' ' + entries[10] + ' ' + groName
            out1 += ' ' + entries[1] + ' ' + str(new_index) + ' ' + str(
                length) + ' ' + str(weight) + ' 0' + entries[5] + ' 0' + entries[6] + "\n"
        else:
            print(pdbFile, "does not exist! Go figure...")

    return out, out1, Missing_pdb, Missing_count


def create_fragment_memories(database, fasta_file, memories_per_position, brain_damage, fragment_length, 
         pdb_dir, index_dir, frag_lib_dir, failed_pdb_list_file, pdb_seqres,
         weight, evalue_threshold, cutoff_identical):
        # set up directories
    
    download_pdb_seqres(pdb_seqres)

    pdb_dir.mkdir(exist_ok=True)
    index_dir.mkdir(exist_ok=True)
    frag_lib_dir.mkdir(exist_ok=True)

    if not pdb_dir.exists() or not frag_lib_dir.exists() or not index_dir.exists():
        print("Can't create necessary directories")
        exit()

    failed_download_pdb_list = []
    if failed_pdb_list_file.exists():
        with failed_pdb_list_file.open() as f:
            failed_download_pdb_list = [line.strip() for line in f]

    #Convert paths to strings
    database = str(database)
    pdb_dir=str(pdb_dir)+'/'
    index_dir=str(index_dir)+'/'
    frag_lib_dir=str(frag_lib_dir)+'/'
    failed_pdb_list_file=str(failed_pdb_list_file)
    pdb_seqres=str(pdb_seqres)

    handle = open(fasta_file, "r")
    N_mem=memories_per_position
    residue_base=0


    for record in SeqIO.parse(handle, "fasta"):
        # Part I, BLAST fragments
        if(len(record.seq) < fragment_length):
            print("Exception::query sequence is shorter than " + str(fragment_length) + " residues. Exit.")
            sys.exit()

        query = str(record.name)[0:4]
        print('processing sequence:', record.name)


        # FRAGMENT GENERATION LOOP
        iterations = len(record.seq) - fragment_length + \
            1  # number of sliding windows

        with ThreadPoolExecutor(max_workers=12) as executor:
            futures={}
            for i in range(1, iterations + 1):
                future = executor.submit(process_window, i, record, fragment_length, evalue_threshold, database, residue_base)
                futures[future] = i-1
            results = [None] * len(futures)
            for future in as_completed(futures):
                i = futures[future]
                results[i] = future.result()
            
        # Writing the results to the match file in the order of execution
        with open('prepFrags.match', 'w') as match:
            match.write(query + "\n")
            for result in results:
                for line in result:
                    match.write(line + '\n')


        # loop2 close
        with open('prepFrags.match', 'r') as match:
            # list unique PDB IDs for downloading later
            matchlines = list()
            keys = {}
            for line in match.readlines():
                matchlines.append(line)
                entries = line.split()
                entry = entries[0]
                if entry[:3] == "pdb":
                    # for example 'pdb|4V12|A'
                    pdbfull = str(entry[4:8]) + str(entry[9:])
                else:
                    pdbfull = str(entry)
                keys[pdbfull] = 1
            unique = list(keys.keys())

        # pdbparse=PDBParser(PERMISSIVE=1)

        # Part II, BLAST the whole sequence to find homologs
        print(record.seq)
        fragment = open('fragment.fasta', 'w')
        fragment.write(str(record.seq))
        fragment.close()
        
        #Download PDBs
        failed_pdb = {pdb_id[0:4].lower():1 for pdb_id in unique}
        failed_pdb.update(download_pdbs([pdb_id[0:4].lower() for pdb_id in unique if pdb_id[0:4].lower() not in failed_download_pdb_list], pdb_dir))
        homo = {pdbID:0 for pdbID in failed_pdb if not failed_pdb[pdbID]}
        homo_count = {pdbID:0 for pdbID in failed_pdb if not failed_pdb[pdbID]}
        update_failed_pdb_list(failed_pdb,failed_pdb_list_file)

        # blast the whole sequence to identify homologs Evalue 0.005
        exeline = "psiblast -num_iterations 1 -word_size 3 -evalue 0.005"
        exeline += " -outfmt '6 sseqid slen bitscore score evalue pident' -matrix BLOSUM62 -db " + \
            database + " -query fragment.fasta"
        print("finding homologs")
        print("executing::: " + exeline)
        homoOut = os.popen(exeline).read()
        homoOut = homoOut.splitlines()  # now an array
        for line in homoOut:
            entries = line.split()
            print("homologues: ", entries)
            if len(entries):
                pdbfull = entries[0]
                pdbID = pdbfull[0:4].lower()
                if brain_damage == 2:
                    identity = float(entries[5])
                    # exclude self(>90% identity)
                    if identity <= cutoff_identical:
                        homo[pdbID] = 1
                        homo_count[pdbID] = 0
                if brain_damage == 0.5:
                    identity = float(entries[5])
                    # check identity, add only self (>90% identity) to homo[]
                    if identity > cutoff_identical:
                        homo[pdbID] = 1
                        homo_count[pdbID] = 0
                else:
                    homo[pdbID] = 1
                    homo_count[pdbID] = 0

        # Part III, Write memories
        count = {}
        for i in range(1, iterations + 1):
            count[str(i)] = 0  # count number of mem per fragments
        fastFile = "./tmp.fasta"
        Missing_pdb={}
        Missing_count=0

        with ThreadPoolExecutor(max_workers=12) as executor:
            futures={}
            for iter,line in enumerate(matchlines):
                future = executor.submit(create_index_files,iter, line, N_mem, brain_damage,count, failed_pdb,homo, homo_count, weight, frag_lib_dir, pdb_dir, index_dir, pdb_seqres)
                futures[future] = iter
            results = [None] * len(futures)
            for future in as_completed(futures):
                i = futures[future]
                results[i] = future.result()
            
        # Writing the results to the match file in the order of execution
        with open('frags.mem', 'w') as LAMWmatch, open('log.mem', 'w') as log_match:
            LAMWmatch.write('[Target]' + "\n")
            LAMWmatch.write("query" + "\n\n" + '[Memories]' + "\n")
            for result in results:
                if result:
                    out, out1, Missing_pdb_out, Missing_count_out = result
                    LAMWmatch.write(out)
                    log_match.write(out1)
                    Missing_pdb.update(Missing_pdb_out)
                    Missing_count+=Missing_count_out

        print("HOMOLOGS:::")
        total_homo_count = 0
        for line in homoOut:
            entries = line.split()
            print("sseqid slen bitscore score evalue pident")
            print(entries)
            entry = entries[0]
            if entry[:3] == "pdb":
                # for example 'pdb|4V12|A'
                pdbfull = str(entry[4:8]) + str(entry[9:])
            else:
                pdbfull = str(entry)
            # pdbfull = entries[0]
            pdbID = pdbfull[0:4].lower()
            if brain_damage == 0 or brain_damage == 2:
                total_homo_count += homo_count[pdbID]
                print("Homolog count =", homo_count[pdbID])

        if brain_damage == 0 or brain_damage == 2:
            print("Total homolog count = ", total_homo_count, round(total_homo_count / iterations, 2))

        print("Memories per position that are fewer than expected:")
        for i in count:
            if count[i] < N_mem:
                print(i, count[i])

        print("Number of blasted PDB: ", len(failed_pdb))
        print("Number of failed downloaded PDB: ", sum(failed_pdb.values()))
        print("Number of PDB with Missing atoms: ", len(Missing_pdb))
        print("Discarded fragments with Missing atoms: ", Missing_count)
        residue_base += len(record.seq)
        for line in homoOut:
            entries = line.split()
            if len(entries):
                entry = entries[0]
                if entry[:3] == "pdb":
                    # for example 'pdb|4V12|A'
                    pdbfull = str(entry[4:8]) + str(entry[9:])
                else:
                    pdbfull = str(entry)
                # pdbfull = entries[0]
                pdbID = pdbfull[0:4].lower()
                print(pdbID)

    # loop1 close

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process some inputs.")
    parser.add_argument("database_prefix")
    parser.add_argument("fasta_file", type=Path)
    parser.add_argument("--memories_per_position", type=int, default=20)
    parser.add_argument("--brain_damage_flag", type=float, default=1)
    parser.add_argument("--frag_length", type=int, default = 9)
    parser.add_argument("--pdb_dir", type=Path, default=Path("PDBs"))
    parser.add_argument("--index_dir", type=Path, default=Path("Indices"))
    parser.add_argument("--frag_lib_dir", type=Path, default=Path("Gros"))
    parser.add_argument("--failed_pdb_list_file", type=Path, default=Path("notExistPDBsList"))
    parser.add_argument("--pdb_seqres", type=Path, default=Path("pdb_seqres.txt"))
    parser.add_argument("--weight", type=int, default=1)
    parser.add_argument("--evalue_threshold", type=int, default=10000)
    parser.add_argument("--cutoff_identical", type=int, default=90)

    args = parser.parse_args()

    create_fragment_memories(args.database_prefix, args.fasta_file, args.N_mem, args.brain_damage_flag, 
         args.frag_length, args.pdb_dir, args.index_dir, args.frag_lib_dir, args.failed_pdb_list_file, args.pdb_seqres,
         args.weight, args.evalue_threshold, args.cutoff_identical)
