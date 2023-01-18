#!/projects/ag-fbz/rnkouame/software/miniconda3/bin/python

import sys

"""
Edit variants format in VCF output. add an info named SVORI and change the SVTYPE to DEL, INS, INV, DUP + switch DUP to INS. 
0- Chrom    1- Pos
2- ID       3- REF
4- ALT      5- QUAL
6- FILTER   7- INFO
8- FORMAT   9- SAMPLE

"""
# transform a list with a specified caracter into a dictionnary
def list_to_dict(liste,separator):
    result= {}
    for element in liste:
        element=element.split(separator)
        try:
            result[element[0]] = element[1]
        except IndexError:
            result[element[0]] = ''
    return result

# transform a dict into a list and separate keys and value with a specified caracter
def dict_to_list(dict,separator):
    result=[]
    for element in dict.keys():
        result.append(element + separator + dict[element])
    return result

def main():
    
    file_path = sys.argv[1]
    output_file = sys.argv[2]
    f_new = open(output_file, "w")
    stop = 0
    info_flag = True
    with open(file_path,"r") as f:
        for line in f:
            if line.startswith('##INFO') and info_flag:
                f_new.write('##INFO=<ID=SVORI,Number=1,Type=String,Description="Original complete description of the SV">\n') 
                info_flag = False
                f_new.write(line)
            elif line.startswith('#'):
                f_new.write(line)
            else:
                stop+=1
                record=line.strip().split('\t')
                record[7]=record[7].split(';')
                info_dict = list_to_dict(record[7],'=')
                info_dict['SVORI'] = info_dict['SVTYPE']
                fine = info_dict['SVTYPE'][0:3]
                if fine =="DUP":
                    fine="INS"
                info_dict['SVTYPE'] = fine

                record[7]=';'.join(dict_to_list(info_dict,'='))
                record = '\t'.join(record)
                f_new.write(record + '\n')


    f_new.close()

if __name__ == "__main__":
    main()
