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
 * Converter for neutrons from ILC extraction line text file to stdhep
 * Author: Anne Schuetz
 */

public class NeutronsToStdhep_SplitFile {

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
					+ ">> java -cp bin:lib/* NeutronsToStdhepLCIO -h / --help");
			System.exit(1);
		}

		//Input and Output files:
		boolean inputfile_set = false;
		boolean outputfile_set = false;
		String input_filename = null;
		String output_filename = null;
		int total_number_of_particles = 0;

		//Find the input and output file:
		for (int i = 0; i < args.length; i++){
			if ( args[i].equals("-h") || args[i].equals("--help")) Usage();
			if ( args[i].equals("-i")){
				input_filename = args[i+1];
				inputfile_set = true;
			}
			if ( args[i].equals("-o")){
				output_filename = args[i+1];
				outputfile_set = true;
			}
		}
		//Check if input and output filenames are set by command line input:
		if (!inputfile_set || !outputfile_set){
			System.out.println("You didn't give an inputfile/outputfile. Please try again!\n");
			Usage();
		}
		//Check if input file with given name exists:
		File inputfile = new File(input_filename);
		if (!inputfile.exists()) {
			System.out.println("Input file " + input_filename + " does not exist!");
			System.exit(1);
		}
		//CountLines: counts lines in file, will exit the program if the file is empty:
		total_number_of_particles = CountLines(input_filename);
		
		//Create the followings for reading in text out of input file:
		FileInputStream neutrons_file = new FileInputStream(inputfile);
		BufferedReader neutrons = new BufferedReader(new InputStreamReader(neutrons_file));
	
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
		int nmax = total_number_of_particles;
		boolean distribution = false; 

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
		}
		//Program exited, if values for p_T or Theta are negative:
		if (nmax < 0) {
			System.out.println("Please give positive values for the maximum number of particles!");
			System.exit(1);
		}
		
		List<Integer> list = new ArrayList<Integer>(); 
		Random randomGenerator = new Random();

		if (distribution == true){
						for(int i = 0; i < total_number_of_particles;++i){
										list.add(i);
						}
						int poisson_n = Poisson(nmax);
						for(int i = 0; i < (total_number_of_particles - poisson_n); ++i){
										int size = list.size();
										int j = randomGenerator.nextInt(size);//gives a random integer in the interval [0;size[
										list.remove(j);
						}
						nmax = poisson_n;
						System.out.println("New nmax = " + nmax);
		}
		else {
						list.clear();
		}

		boolean MoreOutputfiles = false;
		
		if (file_format.equals("stdhep")) {
			/* If given number of particles (nmax) is below total number of particles in input file,
			 * the choice is given to create several output files with each a new set of nmax particles
			 */
			if (nmax < total_number_of_particles) {
				MoreOutputfiles = YesNo_MoreOutputfiles();
			}
			//Call function to convert into stdhep:
			ToStdhep(output_name, neutrons, total_number_of_particles, nmax, list, MoreOutputfiles);
			neutrons_file.close();
		}
		else {
			System.out.println("Unknown file format! Please type output.stdhep!");
			System.exit(1);
		}

	}//end main

	/**
	 * @param outputFilename Name of output file
	 * @param neutrons BufferedReader for input neutrons file
	 * @param tot_num Total number of particles in the input file
	 * @param _nmax Maximum number of particles per output file
	 * @param list List of particle number that are randomly chosen
	 * @param More_outputfiles Boolean choice whether several output files shall be created with each _nmax particles
	 * @param pT_cut_low Value for lower limit of p_T
	 * @param pT_cut_high Value for higher limit of p_T
	 * @param Theta_cut_low Value for lower limit of theta
	 * @param Theta_cut_high Value for higher limit of theta
	 */
	public static void ToStdhep(String outputFilename, BufferedReader neutrons, int tot_num, int _nmax, List<Integer> list, boolean More_outputfiles) {
		
		//Double array of values per line in the input file:
		double[] values = new double[7];

		int _n = 0; //Number of particles
		int _m = 0; //Number of particles
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
					"converted from a neutron text file", 10);
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
			
			//Loop over lines in input file:
			String line;
			while ((line = neutrons.readLine()) != null) {	
				if (!list.isEmpty()){//if the list is empty, then all the particles should be used
								System.out.println("Checking if item "+_m+" exists: " +list.contains(_m));
								if(!list.contains(_m)){
												_m++;
												continue;//if it is not empty, forget about the ones that are not in the list
								}
				}
				
				//If more output files are wished, open new output file and create new event:
				if (_n >= _nmax && More_outputfiles){
					
					//The new output filename has an number extension (_eventnum is counted up for every new output file):
					New_outputFilename = outputFilename + "_" + Integer.toString(_eventnum);
					File NEWoutputfile = new File(New_outputFilename+".stdhep");
					if (NEWoutputfile.exists()) {
						System.out.println("\nAs I wanted to create several output files with each "+ _nmax 
							+ " MCParticles, I was about to create an output file with the filename " 
							+ New_outputFilename + ".stdhep.\n"
							+ "But such a file already exists.\n"
							+ "Please pick another output filename or move the existing file!");
						System.exit(1);
					}
					try {
						//Create a new StdhepWriter with the given output filename:
						w = new StdhepWriter(New_outputFilename+".stdhep", "Stdhep events",
								"converted from a neutrons.txt file", 10);
						w.setCompatibilityMode(false);
						w.writeRecord(new StdhepBeginRun(nevtreq, nevtgen, nevtwrt,
								stdecom, stdxsec, stdseed1, stdseed2));
					}
					//Catch the exception if opening the StdhepWriter for the output file didn't work:
					catch (java.io.IOException e) {
						System.err.println("Error opening file: " + New_outputFilename + ".stdhep");
						e.printStackTrace();
						System.exit(1);
					}
					_n = 0;
					
					/* New arrays for StdhepEvent values for each new output file:
					 * (In this way the arrays always have the right size.)
					 */
					int eff_num = _nmax;
					if (tot_num - (_eventnum-1)*_nmax < _nmax){
						eff_num = tot_num - (_eventnum-1)*_nmax; //The remaining particles
					}
					System.out.println("eff_num = " + eff_num);
					_fst = new int[eff_num];
					_id = new int[eff_num];
					_jmo = new int[2 * eff_num];
					_jda = new int[2 * eff_num];
					_p = new double[5 * eff_num];
					_v = new double[4 * eff_num];
				}
						
				int j = 0;
				StringTokenizer st = new java.util.StringTokenizer(line, " ");
				
				//Store values in the current line into array:
				while (st.hasMoreElements()) {
					values[j++] = Double.valueOf(st.nextToken()).doubleValue();
				}
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
				 
				//If more output files are wished, open new output file and create new event:
				if (_n >= _nmax && More_outputfiles){
					try{
						System.out.println("\n _eventnum = " + _eventnum + ", _n = " + _n + "\n");
						//Write out arrays to the stdhep event:
						StdhepEvent event = new StdhepEvent(_eventnum, _n, _fst, _id,
							_jmo, _jda, _p, _v);
						//Write out stdhep event into the stdhep file:
						w.writeRecord(event);									
						w.writeRecord(new StdhepEndRun(nevtreq, nevtgen, nevtwrt, stdecom,
							stdxsec, stdseed1, stdseed2));
					
						//Write and close the stdhep file:
						System.out.println("\n DONE! Closing file " + New_outputFilename + ".stdhep with " + _n + " MCParticles.");
						w.close();
					}
					//Catch the exception if writing and closing the stdhep file didn't work:
					catch(java.io.IOException ex){
						System.err.println("Error closing file: " + New_outputFilename + ".stdhep");
		                ex.printStackTrace();
		                System.exit(1);
					}
					//_eventnum is counted up for every new output file:
					_eventnum++;
				}
				//Break if maximum number of particles is reached and no more output files are wished:
				if (_n >=_nmax && !More_outputfiles) break;
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
	 * @return Returns boolean for choice, whether more output files should be created or not
	 */
	private static boolean YesNo_MoreOutputfiles() {
		boolean more_outputfiles = false;
		Scanner keyboard = null;

		try {
			System.out
					.println("\n The maximum number of events you have given is smaller than the number of events in the input file. \n"
							+ "Do you want me to create several output files (with each nmax events) until the end of the input file is reached? \n"
							+ "(y/n) \n");
			//Keyboard input is stored in String yes_no:
			keyboard = new Scanner(System.in);
			String yes_no = keyboard.nextLine();

			//Change value of boolean more_outputfiles depending on keyboard input:
			if (yes_no.equals("n")) {
				System.out
						.println("\n Okay, I just convert nmax events and don't care about the rest.");
			} 
			else if (yes_no.equals("y")){
				more_outputfiles = true;
				System.out
						.println("\n Okay, I will create several output files to convert all events of the input file.");
			}
			else {
				System.out.println("\n Please, try again.\n");
				YesNo_MoreOutputfiles();
			}
			return more_outputfiles;
		} finally {
			//End use of keyboard:
			if (keyboard != null) keyboard.close();
		}

	}

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
		System.out.print(String.format("%65s", " ").replace(' ', '*'));
		System.out.printf("%n%-4s%20s%17s%20s%4s%n","****"," ","NeutronsToStdhepLCIO"," ","****");
		System.out.printf("%-4s%57s%4s%n","****"," Converting neutrons.txt files to Stdhep or LCIO ","****");
		System.out.printf("%-4s%18s%21s%18s%4s%n","****"," ","Author:  Anne Schuetz"," ","****");
		System.out.print(String.format("%65s", " ").replace(' ', '*'));
		System.out.println("\n");
	}
	
	private static void Usage() {
		System.out.println("\nNeutronsToStdhepLCIO: \n"
			+ "Application to convert neutrons.txt files to stdhep format or slcio format.\n"
			+ "\nCuts on pT and Theta can be applied in the following way: pTcut_low < pT [GeV] < pTcut_high, and thetacut_low < theta [degrees] < thetacut_high. \n"
			+ "With giving an integer number, the number of particles that are to be converted can be defined.\n"
			+ "By setting the flag -d to true (-d true), the particle selection can be randomized, and the total number of particles can be poisson distributed around the given max. number.\n"
			+ "Passing a run number as an argument will allow to distinguish single simulation files after merging.\n"
			+ "If no run number will be passed the default value 1 will be set as the run number for this file.");
		System.out.println("\nUSAGE: \n"
			+ ">> java -cp bin:lib/*  NeutronsToStdhepLCIO -i PATH/TO/input.txt -o output<.stdhep / .slcio> <more options> \n");
		System.out.println("\nRequired Arguments:\n");
		System.out.printf("%-25s%s%n","-i:","<input.txt file>");
		System.out.printf("%-25s%s%n","-o:","<output filename.stdhep / .slcio>");
		System.out.println("\nOPTIONS:\n");
		System.out.printf("%-25s%s%n","-h / --help:","Usage");
		System.out.printf("%-25s%s%n","-r:","<run number>");
		System.out.printf("%-25s%s%n","-n:","<maximum number of particles that are to be converted>");
		System.out.printf("%-25s%s%n","-d:","<true/false> Boolean to say if the maximum number of particles should be poisson distributed -> random selection of particles");
		System.out.printf("%-25s%s%n","-pl / --ptcut_low:","<lower limit for pT in GeV>");
		System.out.printf("%-25s%s%n","-ph / --ptcut_high:","<higher limit for pT in GeV>");
		System.out.printf("%-25s%s%n","-tl / --thetacut_low:","<lower limit for theta in degree>");
		System.out.printf("%-25s%s%n","-th / --thetacut_high:","<higher limit for theta in degree>");
		System.out.println("\n For example: \n"
			+ ">> java -cp bin:lib/* NeutronsToStdhepLCIO -i neutrons.txt -o neutrons.slcio -r 2 -n 3000 -pl 0.01 -ph 1 -tl 0.2 -th 30");
		System.exit(0);
	}//end Usage()

	private static int Poisson(int mean){
					double L = Math.exp(-mean);
					double p = 1.0;
					int k = 0;

					do {
									k++;
									p *= Math.random();
					} while (p > L);

					return k - 1;
	}

}//end NeutronsToSthepLCIO class


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
		pos[0] = qualities[0];//x (mm)
		pos[1] = qualities[1];//y (mm)
		pos[2] = 9300.0;//z (mm) (neutrons at IP)
		cos = new double[3];
		cos[0] = qualities[2];//directional cosine between x direction and z-axis
		cos[1] = qualities[3];//directional cosine between y direction and z-axis
		cos[2] = qualities[4];//directional cosine between z direction and z-axis
		mom_tot = qualities[5];
		time = qualities[6]; //time after neutrons created 

		//Default values for particle mass, pdg, charge and state of simulation:
		GeneratorStatus = 1; //These particles are all generated.
		pdg = 2112; //PDG Monte Carlo numbering scheme
		mass = 0.9395654133; //neutron mass in 
		charge = 0;
		
	}
	/* Initialise the particle by calculating the momentum, p_T and theta.
	 * Additionally, check whether particle is an electron or positron.
	 */
	public void initialise(){
		energy = calculateEnergy(mom_tot);

		mom = new double[3];
		mom = calculateMom(mom_tot,cos);

		boolean positronline = false;//CHANGE HERE WETHER IT IS THE POSITRON LINE OR NOT
		if (positronline) {
			ChangeToPositronBeam();	
		}
	}
	private double[] calculateMom(double mom_tot, double[] cos){
		double[] p = new double[3]; 
		p[0] = mom_tot * cos[0];
		p[1] = mom_tot * cos[1];
		p[2] = mom_tot * cos[2];
		return p;
	}
	private double calculateEnergy(double mom_tot){
		double E = Math.sqrt( Math.pow(mom_tot,2) + Math.pow(mass,2) );
		return E;
	}
	//Switch sign of x coordinate and momenta
  public void ChangeToPositronBeam(){
		pos[0] *= -1.0;
		pos[2] *= -1.0;
		cos[0] *= -1.0;
		mom = calculateMom(mom_tot,cos);
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
	private double[] cos;
	//private double[] beta;
	private double[] pos;
	private double energy;
	private int pdg = 2112; //PDG number of an neutron
	private int charge; 
	private double mass = 0.9395654133; //neutron mass
	private double time;
	private int GeneratorStatus;
}
