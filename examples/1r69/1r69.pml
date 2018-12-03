sequence = "SISSRVKSKRIQLGLNQAELAQKVGTTQQSIEQLENGKTKRPRFLPELASALGVSVDWLLNGT"
for aa in sequence: cmd._alt(str.lower(aa))
alter (all), resi=str(int(resi)-1)
save 1r69.pdb
quit
