import hep.io.stdhep.StdhepBeginRun;
import hep.io.stdhep.StdhepEndRun;
import hep.io.stdhep.StdhepEvent;
import hep.io.stdhep.StdhepWriter;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.*;
//import java.util.Scanner;
//import java.util.StringTokenizer;

//import java.util.Random;

/**
 * Converter for muons from spoilers text file to stdhep
 * Author: Anne Schuetz
 */

public class MuonsToStdhep {

	/**
	 * @param args String array of command line input
	 * @throws IOException
	 */
	public static void main(String[] args) throws IOException {
		Start();
		if (args.length < 1) Usage();
		if (args.length %2 != 0 && !(args[0].equals("-h") || args[0].equals("--help"))){
			System.out.println("Please check your arguments!\n"
					+ "I guess you forgot to set a flag... Type for the USAGE:\n"
					+ ">> java -cp bin:lib/* MuonsToStdhepLCIO -h / --help");
			System.exit(1);
		}

		//Input and Output files:
		boolean inputfiles_set = false;
		boolean outputfile_set = false;
		List<String> input_filenames = new ArrayList<String>(); 
		String output_filename = null;
		List<Integer> number_of_particles = new ArrayList<Integer>(); 

		//Find the input and output file:
		for (int i = 0; i < args.length; i++){
			if ( args[i].equals("-h") || args[i].equals("--help")) Usage();
			if ( args[i].equals("-i")){
				int j = 1;
				while( args[i+j] != null
					&& !args[i+j].equals("-h") && !args[i+j].equals("--help")
					&& !args[i+j].equals("-o")
					&& !args[i+j].equals("-r")
					&& !args[i+j].equals("-n")
					&& !args[i+j].equals("-d")
					&& !args[i+j].equals("-weights")
 					){
					input_filenames.add( args[i+j] );
					j++;
				}
				inputfiles_set = true;
			}
			if ( args[i].equals("-o")){
				output_filename = args[i+1];
				outputfile_set = true;
			}
		}
		//Check if input and output filenames are set by command line input:
		if (!inputfiles_set || !outputfile_set){
			System.out.println("You didn't give all inputfiles/outputfile. Please try again!\n");
			Usage();
		}
		//CountLines: counts lines in file, will exit the program if the file is empty:
		//Check if input file with given name exists:
		List< File > inputfiles = new ArrayList< File >(); 
		for(int file = 0; file < input_filenames.size(); file++){
			File new_file =  new File(input_filenames.get(file));
    			//System.out.println(input_filenames.get(file));
			if (!new_file.exists()) {
				System.out.println("Input file " + input_filenames.get(file) + " does not exist!");
				System.exit(1);
			}
			BufferedReader br = new BufferedReader(new FileReader( input_filenames.get(file) ));     
			if (br.readLine() == null) {
    				System.out.println(input_filenames.get(file) + " is empty");
				continue;
			}
			else{
				inputfiles.add( new_file );
				number_of_particles.add( CountLines(input_filenames.get(file)) );
			}
		}
    		
		System.out.println("File contents:");
		for(int file = 0; file < inputfiles.size(); file++){
    			System.out.println(inputfiles.get(file).getName()+": "+number_of_particles.get(file));
		}
		//System.exit(1);
		
		//Create the followings for reading in text out of input file:
		List< FileInputStream > muons_files = new ArrayList< FileInputStream >();
		List< BufferedReader > muons = new ArrayList< BufferedReader >();
		for(int file = 0; file < inputfiles.size(); file++){
			muons_files.add( new FileInputStream(inputfiles.get(file)) );
			muons.add( new BufferedReader(new InputStreamReader(muons_files.get(file))) );
		}
	
		//Check if output file already exists -> if yes, program is exited:
		File outputfile = new File(output_filename);
		if (outputfile.exists()) {
			System.out.println("Output file " + output_filename + " already exists!\n"
					+ "Please pick another output filename!");
			System.exit(1);
		}
		//Split output filename into name and format ending:
		int dot = output_filename.lastIndexOf(".");
		String file_format = output_filename.substring(dot+1);
		String output_name = output_filename.substring(0,dot);

		//Default values for the optional arguments:	
		int runnum = 1;
		int nmax = 0;
		boolean distribution = false; 
		List< Double > weights = new ArrayList< Double >(); 
		//double weight_ABE  = 0.0;
		//double weight_PC7  = 0.0;
		//double weight_PC6  = 0.0;
		//double weight_AB5  = 0.0;
		//double weight_PC5A = 0.0;
		//double weight_PC5  = 0.0;
		//double weight_SP4  = 0.0;
		//double weight_PC2  = 0.0;
		//double weight_AB3  = 0.0;
		//double weight_PC1  = 0.0;
		//double weight_SP2  = 0.0;

		//Save values for the optional input values:
		for (int i = 0; i < args.length; i++){
			if ( args[i].equals("-r")){
				runnum = Integer.parseInt(args[i+1]);
			}
			if ( args[i].equals("-n")){
				nmax = Integer.parseInt(args[i+1]);
			}
			if ( args[i].equals("-d")){
				distribution = Boolean.parseBoolean(args[i+1]);
			}
			if ( args[i].equals("-weights")){
				int j = 1;
				while( args[i+j] != null
					&& !args[i+j].equals("-h") && !args[i+j].equals("--help")
					&& !args[i+j].equals("-o")
					&& !args[i+j].equals("-r")
					&& !args[i+j].equals("-n")
					&& !args[i+j].equals("-d")
					&& !args[i+j].equals("-i")
 					){
					weights.add( Double.parseDouble(args[i+j]) );
    					//System.out.println(weights.get( weights.size()-1 ));
					if(j == input_filenames.size()) break;
					j++;
				}
			}
		}
		
		List< List<Integer> > lists = new ArrayList< List<Integer> >(); 
		Random randomGenerator = new Random();
		int new_nmax = 0;

		if (distribution == true){
			for(int file = 0; file < inputfiles.size(); file++){
				List<Integer> list = new ArrayList< Integer >(); 
				lists.add( list );
				for(int i = 0; i < number_of_particles.get(file);++i){
					lists.get(file).add(i);
				}
				int random_n = (int)( (randomGenerator.nextGaussian()*0.08 + 1)* nmax*weights.get(file) );
				if(file==0) System.out.println("Number of muons per file:");
				System.out.println(random_n);
				if(number_of_particles.get(file) - random_n < 0 ){
					System.out.println("Not enough muons in file: " + inputfiles.get(file).getName());
					System.exit(1);
				}
				for(int i = 0; i < (number_of_particles.get(file) - random_n); ++i){
					int size = lists.get(file).size();
					int j = randomGenerator.nextInt(size);//gives a random integer in the interval [0;size[
					lists.get(file).remove(j);
				}//The lists contain in the end numbers of those particles that are to be taken
			new_nmax += random_n;
			}
			System.out.println("New nmax = " + new_nmax);
		}
		else {
			lists.clear();
		}
		//System.exit(1);
		
		if (file_format.equals("stdhep")) {
			/* If given number of particles (nmax) is below total number of particles in input file,
			 * the choice is given to create several output files with each a new set of nmax particles
			 */
			//Call function to convert into stdhep:
			ToStdhep(output_name, muons, number_of_particles, new_nmax, lists);
			for(int file = 0; file < inputfiles.size(); file++){
				muons_files.get(file).close();
			}
		}
		else {
			System.out.println("Unknown file format! Please type output.slcio or output.stdhep!");
			System.exit(1);
		}

	}//end main

	/**
	 * @param outputFilename Name of output file
	 * @param muons BufferedReader for input muons file
	 * @param tot_num Total number of particles in the input file
	 * @param _nmax Maximum number of particles per output file
	 * @param list List of particle number that are randomly chosen
	 * @param More_outputfiles Boolean choice whether several output files shall be created with each _nmax particles
	 * @param pT_cut_low Value for lower limit of p_T
	 * @param pT_cut_high Value for higher limit of p_T
	 * @param Theta_cut_low Value for lower limit of theta
	 * @param Theta_cut_high Value for higher limit of theta
	 */
	public static void ToStdhep(String outputFilename, List< BufferedReader > muons, List<Integer> tot_num, int _nmax, List< List<Integer> > lists) {
		
		//Double array of values per line in the input file:
		double[] values = new double[8];

		int _n = 0; //Number of particles
		int _eventnum = 1; //For counting up output files

		
		//Arrays for StdhepEvent values:
		int[] _fst = new int[_nmax];
		int[] _id = new int[_nmax];
		int[] _jmo = new int[2 * _nmax];
		int[] _jda = new int[2 * _nmax];
		double[] _p = new double[5 * _nmax];
		double[] _v = new double[4 * _nmax];
				
		//Dummy values:
		int nevtreq = 1;
		int nevtgen = 1;
		int nevtwrt = 1;
		float stdecom = 2.F;
		float stdxsec = 2.F;
		double stdseed1 = 137.;
		double stdseed2 = 138.;

		StdhepWriter w = null;
		String New_outputFilename = outputFilename;
		
		try {
			//Create a new StdhepWriter with the given output filename:
			w = new StdhepWriter(outputFilename+".stdhep", "Stdhep events",
					"converted from a muon text file", 10);
			w.setCompatibilityMode(false);
			w.writeRecord(new StdhepBeginRun(nevtreq, nevtgen, nevtwrt,
					stdecom, stdxsec, stdseed1, stdseed2));
		}	
		//Catch the exception if opening the output file didn't work:
		catch (java.io.IOException e) {
			System.err.println("Error opening file: " + outputFilename + ".stdhep");
			e.printStackTrace();
			System.exit(1);
		}

		try {
			_n = 0;
			_eventnum = 1;
			
			for(int file = 0; file < tot_num.size(); file++){
				//Loop over lines in input file:
				String line;
				int _m = 0; //Number of particles
				while ((line = muons.get(file).readLine()) != null) {	
					if (!lists.get(file).isEmpty()){
						//System.out.println("Checking if item "+_m+" exists: " +lists.get(file).contains(_m));
						if(!lists.get(file).contains(_m)){
							_m++;
							continue;//if it is not empty, forget about the ones that are not in the list
						}
					}
					else{//if the list is empty, skip this file
						break;
					}
						
					int j = 0;
					StringTokenizer st = new java.util.StringTokenizer(line, " ");
				
					//Store values in the current line into array:
					while (st.hasMoreElements()) {
						values[j++] = Double.valueOf(st.nextToken()).doubleValue();
					}
					//for(int v = 0; v < values.length; v++){
					//	System.out.println( values[v] );
					//}
					//System.exit(1);
					
					//Initialise a MCParticle with these values:
					Particle p = new Particle(values);
					p.initialise();
				
				
					//Store the qualities of the particle into arrays for the stdhep file:
					_fst[_n] = 1; // final state particle
					_id[_n] = (int) p.getPDG();
					_p[0 + 5 * _n] = p.getMomentum()[0]; // px
					_p[1 + 5 * _n] = p.getMomentum()[1]; // py
					_p[2 + 5 * _n] = p.getMomentum()[2]; // pz
					_p[3 + 5 * _n] = p.getEnergy(); // E
					_p[4 + 5 * _n] = p.getMass(); // mass
					_v[0 + 4 * _n] = p.getPosition()[0]; // x
					_v[1 + 4 * _n] = p.getPosition()[1]; // y
					_v[2 + 4 * _n] = p.getPosition()[2]; // z
					_v[3 + 4 * _n] = p.getTime(); // creation time
				
					//Increment the number of particles in this event
					_n++;
					_m++;
					//System.out.println(_n);
					//System.out.println(_m);
				} 
			}
			
			try{
				System.out.println("\n _eventnum = " + _eventnum + ", _n = " + _n + "\n");
				//Write out the (last) arrays to the (last) stdhep event:
				StdhepEvent event = new StdhepEvent(_eventnum, _n, _fst, _id,
						_jmo, _jda, _p, _v);
				//Write out the (last) stdhep event into the (last) stdhep file:
				w.writeRecord(event);
				w.writeRecord(new StdhepEndRun(nevtreq, nevtgen, nevtwrt, stdecom,
					stdxsec, stdseed1, stdseed2));
				
				//Write and close the (last) stdhep file:
				System.out.println("\n DONE! Closing file " + New_outputFilename + ".stdhep with " + _n + " MCParticles.");
				w.close();
			}
			//Catch the exception if writing and closing the stdhep file didn't work:
			catch(java.io.IOException ex){
				System.err.println("Error opening/closing file: " + New_outputFilename + ".stdhep");
                		ex.printStackTrace();
                		System.exit(1);
			}
			
		} catch (java.io.IOException ex) {
			ex.printStackTrace();
		}

	}//end ToStdhep()

	/**
	 * @param filename Name of file for which lines are to be counted
	 * @return Returns int number of lines
	 * @throws IOException
	 */
	private static int CountLines(String filename) throws IOException {
		//Initialise reading file content:
		BufferedReader reader = new BufferedReader(new FileReader(filename));
		boolean empty = true;
		try {
			int lines = 0;
			//Count up lines till the end of the file:
			while (reader.readLine() != null) {
				empty = false;
				lines++;
			}
			//If file not empty, return the number of lines:
			return (lines == 0 && !empty) ? 1 : lines;
		} finally {
			if (empty) {
				System.out.println("Input file is empty!");
				System.exit(1);
			}
			reader.close();
		}
	}

	private static void Start() {
		System.out.print(String.format("%46s", " ").replace(' ', '*'));
		System.out.printf("%n%-4s%12s%13s%13s%4s%n","****"," ","MuonsToStdhep"," ","****");
		System.out.printf("%-4s%38s%4s%n","****"," Converting muons.txt files to Stdhep ","****");
		System.out.printf("%-4s%8s%21s%9s%4s%n","****"," ","Author:  Anne Schuetz"," ","****");
		System.out.print(String.format("%46s", " ").replace(' ', '*'));
		System.out.println("\n");
	}
	
	private static void Usage() {
		System.out.println("\nMuonsToStdhep: \n"
			+ "Application to convert muons.txt files to stdhep format.\n"
			+ "With giving an integer number, the number of particles that are to be converted can be defined.\n"
			+ "By setting the flag -d to true (-d true), the particle selection can be randomized, and the total number of particles can be Gaussian distributed around the given max. number.\n"
			+ "Passing a run number as an argument will allow to distinguish single simulation files after merging.\n"
			+ "If no run number will be passed the default value 1 will be set as the run number for this file.");
		System.out.println("\nUSAGE: \n"
			+ ">> java -cp bin:lib/*  MuonsToStdhep -i PATH/TO/input.txt -o output<.stdhep> <more options> \n");
		System.out.println("\nRequired Arguments:\n");
		System.out.printf("%-25s%s%n","-i:","<List of input.txt files>");
		System.out.printf("%-25s%s%n","-o:","<output filename.stdhep / .slcio>");
		System.out.println("\nOPTIONS:\n");
		System.out.printf("%-25s%s%n","-h / --help:","Usage");
		System.out.printf("%-25s%s%n","-r:","<run number>");
		System.out.printf("%-25s%s%n","-n:","<maximum number of particles that are to be converted>");
		System.out.printf("%-25s%s%n","-d:","<true/false> Boolean to say if the maximum number of particles should be Gaussian distributed -> random selection of particles");
		System.out.printf("%-25s%s%n","-weights:","<doubles> Weights of different muon inputfiles");
		System.out.println("\n For example: \n"
			+ ">> java -cp bin:lib/* MuonsToStdhep -i muons1.txt muons2.txt -o muons.stdhep -r 2 -n 3000 -weights 0.02 0.23");
		System.exit(0);
	}//end Usage()

	//private static int Poisson(int mean){
	//	double L = Math.exp(-mean);
	//	double p = 1.0;
	//	int k = 0;

	//	do {
	//		k++;
	//		p *= Math.random();
	//	} while (p > L);

	//	return k - 1;
	//}

}//end MuonsToSthepLCIO class


class Particle{
	/**
	 * @param qualities Array of values for the MCParticle energy, beta and vertex
	 */
	public Particle(double[] qualities){
		//Take array of particle qualities and store them in class member vectors:

		//energy = qualities[0];
		//beta = new double[3];
		//beta[0] = qualities[1];
		//beta[1] = qualities[2];
		//beta[2] = qualities[3];
		pos = new double[3];
		pos[0] = qualities[0]*10.0;//x (mm)
		pos[1] = qualities[1]*10.0;//y (mm)
		pos[2] = -10000.0;//z (mm) (muons at IP)
		angle = new double[3];
		angle[0] = qualities[2];//angle between x direction and z-axis
		angle[1] = qualities[3];//angle between y direction and z-axis
		angle[2] = 0.0;//angle between z direction and z-axis
		mom_tot = qualities[4];
		time = -qualities[5]; //time difference between muons and beam (negative sign because negative time in text file means later than beam)
		chargesign= qualities[6]; //charge sign

		//Default values for particle mass, pdg, charge and state of simulation:
		GeneratorStatus = 1; //These particles are all generated.
		pdg = 13; //PDG Monte Carlo numbering scheme
		mass = 0.1056583745; //muon mass in 
		charge = -1;
		
	}
	/* Initialise the particle by calculating the momentum, p_T and theta.
	 * Additionally, check whether particle is an electron or positron.
	 */
	public void initialise(){
		energy = calculateEnergy(mom_tot);

		mom = new double[3];
		mom = calculateMom(mom_tot,angle);

		boolean positronline = true;//CHANGE HERE WETHER IT IS THE POSITRON LINE OR NOT
		if (positronline) {
			ChangeToPositronBeam();	
		}
		if (chargesign > 0) {
			ChangeToPosMuon();
		}
	}
	private double[] calculateMom(double mom_tot, double[] angle){
		double[] p = new double[3]; 
		p[0] = mom_tot * Math.sin(angle[0]);
		p[1] = mom_tot * Math.sin(angle[1]);
		p[2] = Math.sqrt( Math.pow(mom_tot,2) - ( Math.pow(p[0],2) + Math.pow(p[1],2) ) );
		return p;
	}
	private double calculateEnergy(double mom_tot){
		double E = Math.sqrt( Math.pow(mom_tot,2) + Math.pow(mass,2) );
		return E;
	}
	//Switch PDG number, energy value to correspondent values
	public void ChangeToPosMuon(){
		pdg *= -1;
		energy *= -1.0D;
	}
	//Switch sign of x coordinate and angles
	public void ChangeToPositronBeam(){
		pos[0] *= -1.0;
		pos[2] *= -1.0;
		angle[0] *= -1.0;
		mom = calculateMom(mom_tot,angle);
		mom[2] *= -1.0;
	}
	public final double[] getMomentum(){
		return mom;
	}
	public final double[] getPosition(){
		return pos;
	}
	public final double getEnergy(){
		return energy;
	}
	public final int getPDG(){
		return pdg;
	}
	public final int getCharge(){
		return charge;
	}
	public final double getMass(){
		return mass;
	}
	public final double getTime(){
		return time;
	}
	public final int getGeneratorStatus(){
		return GeneratorStatus;
	}
	
	private double mom_tot;
	private double[] mom;
	private double[] angle;
	//private double[] beta;
	private double[] pos;
	private double energy;
	private int pdg = 13; //PDG number of an muon
	private double chargesign; //muon charge
	private int charge; 
	private double mass = 0.1056583745; //muon mass
	private double time;
	private int GeneratorStatus;
}
