sequence = "SISSRVKSKRIQLGLNQAELAQKVGTTQQSIEQLENGKTKRPRFLPELASALGVSVDWLLNGT"
for aa in sequence: cmd._alt(str.lower(aa))
alter (all), resi=str(int(resi)-1)
save extended.pdb
quit
