mol load pdb movie.pdb


graphics 0 line {-1000 15 15} {1000 15 15} width 2
graphics 0 line {-1000 15 -15} {1000 15 -15} width 2
graphics 0 line {-1000 -15 15} {1000 -15 15} width 2
graphics 0 line {-1000 -15 -15} {1000 -15 -15} width 2


while {[molinfo top get numreps] > 0} {mol delrep 0 top}

mol representation NewCartoon
mol material AOChalky
#mol modmaterial 0 0 AOChalky


while {[molinfo top get numreps] > 0} {mol delrep 0 top}
mol color ColorID 1
mol selection resid 1 to 90
mol addrep top

mol color ColorID 2
mol selection resid 91 to 115
mol addrep top


mol color ColorID 3
mol selection resid 116 to 146
mol addrep top


mol color ColorID 4
mol selection resid 147 to 169
mol addrep top

mol color ColorID 7
mol selection resid 170 to 196
mol addrep top

mol color ColorID 10
mol selection resid 197 to 223
mol addrep top

mol color ColorID 0
mol selection resid 224 to 248
mol addrep top

mol color ColorID 11
mol selection resid 248 to 271
mol addrep top

axes location off
display projection orthographic
display cuedensity 0
color Display Background white

user add key q {rotate x by 90}
user add key w {rotate x by -90}

puts "ahora anda lindo todo!!!"

animate style Once
animate goto 0
display resetview

rotate x by 90.000000
rotate x by 90.000000
rotate x by 90.000000

# translate by 0.0 1.0 0.0
# scale by 0.833000
# scale by 0.9
animate speed 1.0

mol smoothrep 0 0 10
mol smoothrep 0 1 10
mol smoothrep 0 2 10
mol smoothrep 0 3 10
mol smoothrep 0 4 10
mol smoothrep 0 5 10
mol smoothrep 0 6 10
mol smoothrep 0 7 10
mol smoothrep 0 8 10

#animate forward
#animate goto 450
#source "~/Downloads/take_picture.tcl"
#take_picture
#animate goto 200
#render Tachyon frame200.dat "'/Applications/VMD 1.9.2.app/Contents/vmd/tachyon_MACOSXX86'" -aasamples 12 %s -format TARGA -o frame200.tga -res 2000 2000
#animate goto 450
#render Tachyon frame450.dat "'/Applications/VMD 1.9.2.app/Contents/vmd/tachyon_MACOSXX86'" -aasamples 12 %s -format TARGA -o frame450.tga -res 2000 2000
#exit
