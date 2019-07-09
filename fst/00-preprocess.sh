#!/bin/bash

cat sample_index.txt | sed 's/ACB\|ASW\|ESN\|MSL\|GWD\|LWK\|YRI/Africa/g' | sed 's/IBS\|GBR\|FIN\|TSI\|CEU/WestEurasia/g' | sed 's/ITU\|STU\|BEB\|PJL\|GIH/SouthAsia/g' | sed 's/KHV\|CDX\|CHS\|JPT\|CHB/EastAsia/g' | sed 's/PEL\|CLM\|PUR\|MXL/America/g' | sed 's/NaN/EastAsia/g' > sample_index_manualassign.txt
