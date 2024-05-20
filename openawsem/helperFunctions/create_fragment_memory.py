import argparse
from pathlib import Path
import sys
import os
from .IndexPdb import *
from .Pdb2GroLib import *
from Bio import SeqIO
import subprocess
import argparse
from pathlib import Path


def create_fragment_memories(database, fasta_file, memories_per_position, brain_damage, fragment_length, 
         pdb_dir, index_dir, frag_lib_dir, failed_pdb_list_file, pdb_seqres,
         weight, evalue_threshold, cutoff_identical):
        # set up directories
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

    # set up out
    LAMWmatch = open('frags.mem', 'w')
    LAMWmatch.write('[Target]' + "\n")
    LAMWmatch.write("query" + "\n\n" + '[Memories]' + "\n")
    log_match = open('log.mem', 'w')

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
        match = open('prepFrags.match', 'w')
        match.write(query + "\n")

        # FRAGMENT GENERATION LOOP
        iterations = len(record.seq) - fragment_length + \
            1  # number of sliding windows
        for i in range(1, iterations + 1):  # loop2, run through all sliding windows in a given chain
            rangeStart = i - 1
            rangeEnd = i + fragment_length - 1
            print("window position:::" + str(i))
            subrange = str(record[rangeStart:rangeEnd].seq)
            print("fragment subrange:::" + subrange)
            fragment = open('fragment.fasta', 'w')
            fragment.write(subrange)
            fragment.close()  # a temporary file for BLAST

            # submit PSI-BLAST, run "psiblast -help" for more details of output
            # format (outfmt)
            exeline = "psiblast -num_iterations 5 -comp_based_stats 0 -word_size 2 -evalue " + \
                str(evalue_threshold)
            #exeline+=" -outfmt '6 sseqid qstart qend sstart send qseq sseq length gaps bitscore evalue' -matrix PAM30 -threshold 9 -window_size 0"
            exeline += " -outfmt '6 sseqid qstart qend sstart send qseq sseq length gaps bitscore evalue' -matrix BLOSUM62 -threshold 9 -window_size 0"
            exeline += " -db " + str(database) + " -query fragment.fasta"
            print("executing:::" + exeline)

            psiblastOut = os.popen(exeline).read()
            psiblastOut = psiblastOut.splitlines()  # now an array
            N_blast = len(psiblastOut)
            # print psiblastOut
            if psiblastOut[N_blast - 1] == 'Search has CONVERGED!':
                N_blast = N_blast - 2  # exclude last two lines for the text

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
            print("Number of searched PDBs:  ", N_blast, N_start_new)

            # convert psiblastOut to a list, sorted by evalue
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
                #print(list_tmp)
            psilist.sort(key=lambda x: x[10])

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
                out += '\n'
                gaps = this[8]
                if(gaps == '0'):
                    match.write(out)  # skip gapped alignments
        # loop2 close
        match.close()
        match = open('prepFrags.match', 'r')  # match is read-only now

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
        homo = {}
        failed_pdb = {}
        homo_count = {}
        for pdbfull in unique:
            pdbID = pdbfull[0:4].lower()
            pdbIDsecond = pdbfull[1:2].lower()
            pdbIDthird = pdbfull[2:3].lower()
            chainID = pdbfull[4:5].lower()
            failed_pdb[pdbID] = 0
            homo[pdbID] = 0
            homo_count[pdbID] = 0

            if pdbID in failed_download_pdb_list:
                failed_pdb[pdbID] = 1
                print(":::Cannot build PDB for PDB ID, skipped:" + pdbID.upper())
                continue
            # download PDBs if not exist    ##from script 'pdbget' (original author
            # unknown)
            if not os.path.isfile(pdb_dir + pdbID.upper() + ".pdb"):
                exeline = "wget ftp://ftp.wwpdb.org/pub/pdb/data/structures/divided/pdb/"
                exeline += pdbIDsecond + pdbIDthird + "/pdb" + pdbID + ".ent.gz"
                os.system(exeline)
                print(exeline)
                os.system("nice gunzip pdb" + pdbID + ".ent.gz; mv pdb" +
                        pdbID + ".ent " + pdb_dir + pdbID.upper() + ".pdb")
            if not os.path.isfile(pdb_dir + pdbID.upper() + ".pdb"):
                failed_pdb[pdbID] = 1
                print(":::Cannot build PDB for PDB ID, failed to download:" + pdbID.upper())
                os.system(f"echo '{pdbID}' >> {failed_pdb_list_file}")

            # exit()

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
        iter = 0
        count = {}
        for i in range(1, iterations + 1):
            count[str(i)] = 0  # count number of mem per fragments
        Missing_count = 0
        Missing_pdb = {}
        fastFile = "./tmp.fasta"

        for line in matchlines:  # loop3
            iter += 1
            if not(iter == 1):
                entries = line.split()
                windows_index_str = entries[11]
                if count[windows_index_str] >= N_mem:
                    continue
                # pdbfull = str(entries[0])
                entry = entries[0]
                if entry[:3] == "pdb":
                    # for example 'pdb|4V12|A'
                    pdbfull = str(entry[4:8]) + str(entry[9:])
                else:
                    pdbfull = str(entry)
                pdbID = pdbfull[0:4].lower()
                pdbIDsecond = pdbfull[1:2].lower()
                pdbIDthird = pdbfull[2:3].lower()
                chainID = pdbfull[4:5].lower()
                groFile = frag_lib_dir + pdbID + chainID + ".gro"
                groName = pdbID + chainID + ".gro"
                pdbFile = pdb_dir + pdbID.upper() + ".pdb"
                indexFile = index_dir + pdbID + chainID + ".index"

                if failed_pdb[pdbID]:
                    continue  # failed-downloaded ones are still in matchlines, need to be ignored
                if brain_damage == 2:
                    if homo[pdbID]:
                        homo_count[pdbID] += 1
                        pass
                    else:
                        print(pdbID, "is not a homolog, discard")
                        continue

                if homo[pdbID]:
                    if brain_damage == 0:
                        print(pdbID, " Using  a homolog.")
                        homo_count[pdbID] += 1
                    if brain_damage == 1:
                        print(pdbID, " is a homolog, discard")
                        continue
                residue_list = entries[6]  # sseq

                res_Start = int(entries[3])
                res_End = int(entries[4])
                print(pdbFile, "start: ", res_Start, "end: ", res_End)
                # Do I have the index file?  If No, write it

                if not os.path.isfile(indexFile):
                    # generate fasta file
                    if not os.path.isfile(pdb_seqres):
                        import urllib
                        print("Need to download pdb_seqres.txt from PDB!")
                        print("Downloading pdb_seqres.txt from ftp://ftp.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt...")
                        url = "ftp://ftp.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt"
                        try:
                            urllib.request.urlretrieve(url, pdb_seqres)
                            print(f"Download complete. Saved to {pdb_seqres}")
                        except urllib.error.URLError as e:
                            print(f"Error downloading file: {e.reason}")
                        except Exception as e:
                            print(f"An error occurred: {e}")
                        print(f"Download complete. Saved to {pdb_seqres}")
                    fastaFile = pdbID + '_' + chainID
                    exeline = "grep -A1 " + fastaFile + " " + pdb_seqres + " > ./tmp.fasta"
                    print("generating fastaFile: ", fastaFile)
                    # p = os.popen(exeline)
                    subprocess.Popen(exeline, shell=True).wait()
                    # p_status = p.wait()
                    if os.path.getsize('./tmp.fasta') > 0:
                        writeIndexFile(fastFile, pdbFile,
                                    indexFile, chainID)
                        print("Writing indexFile: ", indexFile)
                else:
                    print(indexFile, "exist, no need to create.")

                if not os.path.isfile(indexFile):
                    print("Can't create index file, ignore and go on!")
                    continue

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
                    continue
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
                    continue

                if r_list != residue_list:
                    print("Missing residues: ", pdbID + chainID, residue_list, " incomplete: ", r_list)
                    Missing_pdb[pdbID] = 1
                    Missing_count += 1
                    continue

                if os.path.isfile(pdbFile):
                    if not os.path.isfile(groFile):
                        Pdb2Gro(pdbFile, groFile, chainID)
                        print("converting...... " + pdbFile + " --> " + groFile)
                    else:
                        print("Exist " + groFile)
                    count[windows_index_str] += 1

                    length = res_End - res_Start + 1
                    out = groFile + ' ' + entries[1] + ' '  # queue start
                    out += str(new_index) + ' ' + str(length) + ' ' + \
                        str(weight) + "\n"  # frag_seq start
                    LAMWmatch.write(out)
                    out1 = windows_index_str + ' ' + str(count[windows_index_str])
                    out1 += ' ' + entries[9] + ' ' + entries[10] + ' ' + groName
                    out1 += ' ' + entries[1] + ' ' + str(new_index) + ' ' + str(
                        length) + ' ' + str(weight) + ' 0' + entries[5] + ' 0' + entries[6] + "\n"
                    log_match.write(out1)
                else:
                    print(pdbFile, "does not exist! Go figure...")
            # loop3 ends

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

        print("memories per position that is fewer than expected:")
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
