# Creator: Ilian Torres
# Descripton:
# Takes as input the directory of the Main folder called fastqc.
# In the fastqc folder there exists other folders which contain the data txt called fastqc_data.txt.
# This program takes fastqc_data.txt from each folder in fastqc makes a new folder called Sep_fastqc_data.
# In Sep_fastqc_data it creates txt files for each individual data of fastqc_data.txt.

# imports
import os
import argparse
import zipfile
import fnmatch
import shutil
import glob
# Function Make_file_txt
def Make_file_txt(directory):
    i = 1

    name_of_txt = ["Basic Statistics", "Per base sequence quality", "Per tile sequence quality",
                   "Per sequence quality scores", "Per base sequence content","Per sequence GC content", "Per base N content",
                   "Sequence Length Distribution", "Sequence Duplication Levels", "Overrepresented sequences", "Adapter Content","Kmer Content"]
    newpath = os.path.join(directory, "fastqc_data.txt")
    filepath = os.path.join(directory, "Sep_fastqc_data")

    if not os.path.exists(newpath):
        return
    if not os.path.exists(filepath):
        os.mkdir(filepath, 0755)

    C_txt = open(os.path.join(filepath, name_of_txt[0].replace(' ', '_') + '.txt'), "w")

    with open(newpath) as contents:
        for line in contents:
            if ">>END_MODULE" in line:
                C_txt.close()
                try:
                    line = contents.next()
                except StopIteration:
                    break
                if (line.split('\t')[0][2:] == name_of_txt[i]):
                    C_txt = open(os.path.join(filepath, name_of_txt[i].replace(' ', '_') + '.txt'), "w")
                    i = i + 1
                else:
                    print "Error file " + name_of_txt[i] + " not in fastqc_data.txt in directory: " + directory
                    print "Moving to next folder."
                    break
            else:
                C_txt.write(line)

# Main Function
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("directory", help="Creates different txt files for fastqc_data.txt in each folder.", type=str)
    args = parser.parse_args()
    if not os.path.exists(args.directory):
        print "Path directory does not exist.Try again.Press ctr c to quit."
    else:
        rootPath = args.directory
        pattern = '*.zip'
        for root, dirs, files in os.walk(rootPath):
            for filename in fnmatch.filter(files, pattern):
                zipfile.ZipFile(os.path.join(root, filename)).extractall(
                    os.path.join(root, os.path.splitext(filename)[0]))
                os.remove(os.path.join(root,filename))
        for folder in (glob.glob(os.path.join(args.directory, '*'))):
            if os.path.join(args.directory, folder)==os.path.join(args.directory, "notCurrent"):
                shutil.rmtree(os.path.join(args.directory, "notCurrent"))
            else:
                for in_folder in os.listdir(os.path.join(args.directory,folder)):
                    Make_file_txt(os.path.join(args.directory, folder,in_folder))


if __name__ == "__main__":
    main()

