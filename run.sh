# 02/23/2024 13:50:30 CET
(echo -en "\x6C\x1B\x01"; for chr in {1..100}; do tail -c +4 sim.${chr}.bed; done) > merge.bed
# 02/23/2024 13:50:50 CET
cat sim.{1..100}.bim > merge.bim
