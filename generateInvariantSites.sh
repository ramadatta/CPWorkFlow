#NC_013654.1_Escherichia_coli_SE15_DNA_complete_genome  1      3570185  878249  876444  909108  906384  0
#name                                                    start  end      A        T        C        G        N  others
#NC_013654.1_Escherichia_coli_SE15_DNA,_complete_genome  1      4717338  1162728  1161107  1197740  1195763  0

yes A | head -"$1" | tr -d "\n"
yes T | head -"$2" | tr -d "\n"
yes C | head -"$3" | tr -d "\n"
yes G | head -"$4" | tr -d "\n"

echo ""; 
