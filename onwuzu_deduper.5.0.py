#!/usr/bin/python 


##This script removes PCR duplicates from single-end sam files.
##sam file MUST have been previously sorted by RNAME and POS.

# Usage:
    # onwuzu_deduper.py -f <file.sam>
    # onwuzu_deduper.py --file <file.sam>
    # onwuzu_deduper.py -o <prefix>
    # onwuzu_deduper.py --output <prefix>
    # onwuzu_deduper.py -p <paired>
    # onwuzu_deduper.py --paired <paired>
    # onwuzu_deduper.py -u <umi reference file> 
    # onwuzu_deduper.py --umi <umi reference file>


    
import argparse
def get_arguments():
    parser = argparse.ArgumentParser(description=".....")
    parser.add_argument("-f", "--file", help="file path of sam file", required=True, type=str)
    parser.add_argument("-o", "--output", help="prefix for deduped file", required=True, type=str)
    parser.add_argument("-u", "--umi", help="file containing list of UMIs", required=True, type=str)
    parser.add_argument("-p", "--paired", help="file is paired end", required=False, type=str)
    parser.add_argument("-d", "--dups", help="prefix for file containing removed duplicates", required=True, type=str)
    return parser.parse_args()
     
args = get_arguments()



def bit_checker(FLAG):
	'''Takes the bitwise flag and checks it for strandedness. 
	Assumes read it mapped, otherwise returns None.
	Assumes data is single-end.
	Returns "+" on "-", depending on strand'''	
	if (FLAG & 4) !=4:    #check read is mapped
		return None	      #or raise NameError("Read is unmapped")  
	strand = "+"          #starts with + and changes if bit & 16 = 16
	if (FLAG & 16) == 16: #changes strand if bit 16 is set (16 from flag=reverse complimented)
		strand = "-"	
	return strand
	


    
import re 
def cigarettee(CIGAR, POS):
    '''Function to check and adjust POS for soft-clipping. Adjust for only left S'''
    if re.match('^\d+S', CIGAR):
        splitCIGAR = CIGAR.split("S")
        num_clipped = int(splitCIGAR[0])
        adj_POS = POS - num_clipped
        return adj_POS
    else:
        return POS
    

    
if args.paired:
	raise NameError("onwuzu_deduper.py does not process paired-end files at this time")
        

##adds umi reference list to set so as to remove sam lines with umis not in reference list        
umi_dict = set()
with open(args.umi, "r") as umi_ref:
    for line in umi_ref:
        umi_line = line.strip()
        umi_dict.add(umi_line)
    

##and for the main event...    
unique_four = {}    
cur_RNAME = "" 
lines_removed = 0
with open(args.file, 'r') as samfile, open(args.output+"_deduped.sam", 'w') as deduped_file, open(args.dups, 'w') as pcr_dups:
	for line in samfile:
		if line.startswith("@"):
			deduped_file.write(line)
		else:
			splitline = line.split("\t")
			QNAME = splitline[0]
			splitQNAME = QNAME.split(":")
			umi = splitQNAME[7]
            
			if umi not in umi_dict:
				continue
			else:
				FLAG = int(splitline[1])
				RNAME = str(splitline[2])
				POS = int(splitline[3])   
				CIGAR = str(splitline[5])
                
				if RNAME != cur_RNAME:
					cur_RNAME = RNAME
					unique_four.clear()
                
				adj_POS = cigarettee(CIGAR, POS)
				checked_FLAG = bit_checker(FLAG)
                
				check_four = umi+RNAME+str(adj_POS)+str(checked_FLAG)
				if check_four not in unique_four:
                    
					unique_four[check_four] = 0 
					deduped_file.write(line)
				else:
					lines_removed += 1
					pcr_dups.write(line)
					continue
				
			

print(str(lines_removed)+" PCR duplicates have been removed from "+args.file+".")









