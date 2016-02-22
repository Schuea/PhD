#!/usr/bin/python
import os 
import sys, getopt
import subprocess



def main(argv):
    inputfile = ''
    outputfile = ''
    pT_cut = ''
    Theta_cut = ''
    try:
        opts, args = getopt.getopt(argv,"hi:o:p:t:",["ifile=","ofile=","pTcut=","Thetacut="])
    except getopt.GetoptError:
        print 'Usage:'
        print 'test.py -i <inputfile> -o <outputfile> -p <pT cut in GeV> -t <Theta Cut in degrees>'
        sys.exit(2)
    if(len(sys.argv)!=9):
        print 'Usage:'
        print 'test.py -i <inputfile> -o <outputfile> -p <pT cut in GeV> -t <Theta Cut in degrees>'
        sys.exit()
    for opt, arg in opts:
        if opt in ('-h', '--help'):
            print 'Usage:'
            print 'test.py -i <inputfile> -o <outputfile> -p <pT cut in GeV> -t <Theta Cut in degrees>'
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg
        elif opt in ("-o", "--ofile"):
            outputfile = arg
        elif opt in ("-p", "--pTcut"):
            pT_cut = arg
        elif opt in ("-t", "--Thetacut"):
            Theta_cut = arg
    
    if not os.path.isfile(inputfile):
        print("%s does not exist!" % inputfile)
        sys.exit()
    if os.path.isfile(outputfile):
        print("%s already exists! Please pick another output filename!" % outputfile)
        sys.exit()
    print 'Input file is ', inputfile
    print 'Output file is ', outputfile
    print 'pT cut is ', pT_cut
    print 'Theta cut is ', Theta_cut

    output_name = os.path.splitext(outputfile)[0]

    if not(os.path.isfile("pairsTOhepevt/pairToHepEvt.o")):
        process_make = subprocess.Popen(["make"], stdout=subprocess.PIPE, cwd="./pairsTOhepevt")
        print 'Make pairToHepEvt.o ...'
        process_make.wait()

    pairToHepEvt = ("./pairsTOhepevt/pairToHepEvt", inputfile, output_name + ".HEPEvt", pT_cut, Theta_cut)
    process_pairToHepEvt = subprocess.Popen(pairToHepEvt, stdout=subprocess.PIPE)
    print 'Executing pairToHepEvt ...'
    process_pairToHepEvt.wait()

    if not(os.path.isfile("hepevtTOstdhep/HepEVTStdhepGenerator.class")):
        Java_compiler = ('javac', '-classpath', 'hepevtTOstdhep/classpath/*', './hepevtTOstdhep/HepEVTStdhepGenerator.java')
        process_Java_compiler = subprocess.Popen(Java_compiler, stdout=subprocess.PIPE)
        print 'Compiling HepEVTStdhepGenerator.java ...'
        process_Java_compiler.wait()

    text_filename = "ToBeConvertedToStdHep.txt"
    file = open(text_filename,'w')
    file.write(output_name+".HEPEvt")
    file.close()
    print 'Writing HepEvt filename into text file...'

    HepEVTStdhepGenerator = ('java', '-classpath', 'hepevtTOstdhep:hepevtTOstdhep/classpath/*', 'HepEVTStdhepGenerator', text_filename,'.')
    process_HepEVTStdhepGenerator = subprocess.Popen(HepEVTStdhepGenerator, stdout=subprocess.PIPE)
    print 'Executing HepEVTStdhepGenerator...'
    process_HepEVTStdhepGenerator.wait()


    if os.path.isfile(text_filename):
        os.remove(text_filename)
    else:   
        print("Error: %s file not found" % text_filename)

    print 'Conversion done!'

if __name__ == "__main__":
    main(sys.argv[1:])
