for n in {10..1010..50};do ./basic_dgemm $i ; done  > basic_dgemm_test.csv
for n in {10..1010..50};do ./blocked_dgemm $i ; done  > blocked_dgemm_test.csv
for n in {10..1010..50};do ./wrap_dgemm $i ; done  > wrap_dgemm_test.csv
