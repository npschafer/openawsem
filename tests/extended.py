import numpy as np

def create_extended_pdb_from_fasta(filename, output_file_name="output.pdb"):
    # Standard distances and angles for the backbone atoms
    CA_CA =3.8
    CA_C = 1.52/np.sqrt(3/2)
    CA_N = 1.45/np.sqrt(3/2)
    C_O = 1.23
    CA_CB = 1.53/np.sqrt(3/2)
    C_N = 1.33/np.sqrt(3/2)
    t=1/np.sqrt(2)

    # Load the sequence from the FASTA file
    with open(filename, "r") as file:
        lines = file.readlines()
    # Skip the name line and concatenate sequence lines
    sequence = ''.join(line.strip() for line in lines[1:])
    print(sequence)

    # Define coordinates for the first residue
    coord = np.array([[-1*CA_N, 0, -t*CA_N],  #N
                      [0, 0, 0],              #CA
                      [1*CA_C,0,-t*CA_C],     #C
                      [1*CA_C,0,-t*CA_C-C_O], #O
                      [0,-1*CA_CB,t*CA_CB]])  #CB
    coord[:,2]+=t*(CA_C+CA_N)/2

    one_to_three={'A':'ALA', 'C':'CYS', 'D':'ASP', 'E':'GLU', 'F':'PHE',
                  'G':'GLY', 'H':'HIS', 'I':'ILE', 'K':'LYS', 'L':'LEU',
                  'M':'MET', 'N':'ASN', 'P':'PRO', 'Q':'GLN', 'R':'ARG',
                  'S':'SER', 'T':'THR', 'V':'VAL', 'W':'TRP', 'Y':'TYR'}
    with open(output_file_name, "w") as f:
        # For each residue in the sequence, calculate the position of the backbone atoms
        for resid, residue in enumerate(sequence):
            for i,name in enumerate(['N','CA','C','O','CB']):
                pdb_line = f'ATOM  {i*5+1:>5} {name:^4} {one_to_three[residue]:<3} {"A"}{resid:>4}    {coord[i][0]:>8.3f}{coord[i][1]:>8.3f}{coord[i][2]:>8.3f}\n' 
                f.write(pdb_line)
            # f.write("ATOM  %5d  N   %s A %3d    %8.3f%8.3f%8.3f\n" % (i*5+1, residue, i+1, coord[0][0], coord[0][1], coord[0][2]))
            # f.write("ATOM  %5d  CA  %s A %3d    %8.3f%8.3f%8.3f\n" % (i*5+2, residue, i+1, coord[1][0], coord[1][1], coord[1][2]))
            # f.write("ATOM  %5d  C   %s A %3d    %8.3f%8.3f%8.3f\n" % (i*5+3, residue, i+1, coord[2][0], coord[2][1], coord[2][2]))
            # f.write("ATOM  %5d  O   %s A %3d    %8.3f%8.3f%8.3f\n" % (i*5+4, residue, i+1, coord[3][0], coord[3][1], coord[3][2]))
            # f.write("ATOM  %5d  CB  %s A %3d    %8.3f%8.3f%8.3f\n" % (i*5+5, residue, i+1, coord[4][0], coord[4][1], coord[4][2]))
            
            # Calculate new coordinates for the next residue
            coord[:,0]+=CA_CA/np.sqrt(3**2+t**2)*3
            coord[:,1]=-coord[:,1]
            coord[:,2]=-coord[:,2]


        f.write("END\n")

if __name__=="__main__":
    create_extended_pdb_from_fasta("1r69.fasta")