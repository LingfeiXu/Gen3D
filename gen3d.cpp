/*
 * Gen3D, avaliable on sourceforge
 *
 * http://sourceforge.net/projects/gen3d/
 *
 * Any question please send to lxq35@mail.missouri.edu
 */

#include <iostream>
#include<fstream>
#include <iomanip>
#include <math.h>
#include<stdlib.h>
#include <time.h>
#include <string.h>
#include <string>
#include<float.h>

using namespace std;

ofstream mylog;             //log file used to record the scores of models

double dt = 7.0;            //the contact distance threshold
double mymaxd = 20.25;       //the non-contact distance threshold

const int CON = 65;         //con satisfaction threshold
const int NON = 65;			//non-con satisfaction threshold

double maxcon = 0.0;		//record the max contact score
double maxnon = 0.0;		//record the max non-contact score

const int Nt = 10;          //the number of generated vector for each iteration growth step
int addsub = 1;
int flag = addsub;          //used in adaptation
int avgIF = 0;

int step = 50000;           //the step number of adapation
int s = 0;
double Tg = 10.0;		    //temperature for simulated anneanling

int n;						// to determine number of matrices elements
int mcount=0;
double** R;

double rr[Nt][3];           //used in changeDirection, randCord, genetic, growth, adaptation
double P[Nt];               //used in growth, isProbSame
double temp[Nt][3];         //used in growth, calculateCostFunction, adaptation

char* filename = "10";		//normal chromosome 10 is given as input as default
int oz=11;					//given as a default value for oz

int unitBase = 1000000;	    //the base unit in our method have the size of 1MB

// various matrices created
// all except SS are initialized in readInput()
int** S;                    //the contact matrix, used in omitRegions, used in calcScore, calcDistance
int** T;                    //the heatmap matrix which only contains(1,2,3,4) to represents con-satisfied, con-unsatisfied, non-con-satisfied and non-con-unsatisfied, used in omitRegions
double **SS;                //Temperary matrix for ommitted region, used in inMatrix

double** A; //used in genetic
double** B; //used in genetic

double C = 50.0;            //used in calcDistance
double alpha = 1;
double ROOT3 = 1/sqrt(3.0);

int totalIF;



/*
  function to read the input and use data to determine 

  number of elements in matrices and initialize matrices

  no parameters and no return value

  called in beginning of main
*/
void readInput()
{
	ifstream myfile;           // assigned to filename
	string line;               //assigned to line of myfile
	// four variables used to determine min and max value of integers grabbed from file
	unsigned int pos1, pos2;    
	unsigned int min = -1;
	unsigned int max = 0;
	int l, m;                  // assigned to index value of spaces in line


	myfile.open( filename );   // filename is input file
	// while myfile hasn't reached end of file mark
	// used to determine min/max in file to figure out number of elements in matrix
	while( !myfile.eof() )
	{
		getline(myfile, line);

        // assigns index value of first appearance of a space to l
		l = line.find(" ");
		
		// assigns index value of second appearance of a space to m 
		m = line.find(" ", l+1);
		
		// assigns integer value of substring of line, from after l to before m - l
		// file should be looked at to specify what is being grabbed
		pos1 = atoi( line.substr(l+1, m-l-1).c_str() );
		
		// assigns index value to last appearance of a space to l
		l = line.rfind(" ");
		
		// assings integer value of substring of everything in line after l
		pos2 = atoi( line.substr(l+1).c_str() );
		
		// determine the minimum and maximum values based on pos1 and pos2
		if(min>pos1)
			min = pos1;
		if(min>pos2)
			min = pos2;
		if(max<pos1)
			max = pos1;
		if(max<pos2)
			max = pos2;
	}

	// determines number of elements 
	n = (max-min)/unitBase + 1;
	cout <<"n="<< n << endl;
	
	//various matrices are created with n rows
	// S is contact matrix
	S = new int* [n];				//may be we dont need this symmetric array
	T = new int* [n];
	R = new double* [n];
	A = new double* [n];
	B = new double* [n];

	// for some of the matrices make n columns, others use 3 columns
	for(int i=0;i<n;i++)
	{
		S[i] = new int[n];
		T[i] = new int[n];
		R[i] = new double[3];
		A[i] = new double[3];
		B[i] = new double[3];

	}


	// assign all the values of s matrix to 0
	for(int i=0;i<n;i++)
		for(int j=0;j<n;j++)
			S[i][j] = 0;

	// close and reopen file to reset indicater
	myfile.close();
	myfile.open( filename );

	myfile.clear();				//v. v. imp; without this, seekg never works
	myfile.seekg(0, ios::beg);


	// repeat process that was used earlier
	// used to actually get data from file and assign it to matrix S
	while( myfile.good() )
	{
		getline(myfile, line);        // same process as before

		l = line.find(" ");
		m = line.find(" ", l+1);
		pos1 = atoi( line.substr(l+1, m-l-1).c_str() );
		l = line.rfind(" ");
		pos2 = atoi( line.substr(l+1).c_str() );

		// assign data to S matrix
		++S[(pos1-min)/unitBase][(pos2-min)/unitBase];		//IF

	}

	myfile.close();                      // close the file
}




/*	funciton to give input to genetic algorithm, c is 0 or 1 depends on father or mother

	called in main
*/
void printR(int c)
{
	//if c then assign R to A and print R
	//else assign R to B and print R
	if(c)
		for(int i=0;i<n;i++) {
			A[i][0] = R[i][0];A[i][1] = R[i][1];A[i][2]=R[i][2];
			cout << R[i][0] << " " << R[i][1] << " " << R[i][2] << endl;
		}
	else
		for(int i=0;i<n;i++) {
			B[i][0] = R[i][0];B[i][1] = R[i][1];B[i][2]=R[i][2];
			cout << R[i][0] << " " << R[i][1] << " " << R[i][2] << endl;
		}
}



/*  function to clean up various matrices and their data

	no parameters and no return value
	
	called in end of main
*/
void memdel()
{
	for(int i=0;i<n;i++)	// delete values of matrices
	{
		delete []S[i];
		delete []R[i];
		delete []A[i];
		delete []B[i];
	}

	//delete matrices themselves
	delete []S;
	delete []R;
	delete []A;
	delete []B;
}




/*	function to print the contact matrix

	no parameters and no return values

	called in end of main
*/
void printContactMatrix()
{
	//traverse through matrix S, print value and add to avgIF
	for(int i=0;i<n;i++)
	{
		for(int j=0;j<i;j++)
		{
			cout << setw(5) << S[i][j];             // setw(5) used to set width of printed info
			avgIF += S[i][j];
		}
		cout << endl;// << endl;
	}
	totalIF = avgIF;                                // figure total and average IF and print total
	avgIF = avgIF/(n*n);
	cout << "total IF=" << totalIF << endl;
}




/*	function to decide reward or penalized the score base on contact or non-contact satisfaction, used in simualted annealing to calculate energy
	
	parameters int a and return integer 0 or -1

	called in calculateCostFunction with S[i][k] parameters and return
*/
int aG(int a)
{
	return (a==0)?1:-1;                        //if parameter a is zero, return 0 else return -1
}




/*	heaviside function, used to calculated the energy for simulated anneanling

	parameters double n and return integer

	called in calculateCostFunction with dt-rik parameters
*/
int heaviside(double n)
{
	return (n<0)?-1:1;                        //if parameter is less then zero return -1 else return 1
	//return (n<0)?0:1;
}




/*
	function to return distance
	
	parameters includes 6 double values, x1, y1, z1, x2, y2, z2
	
	return is double value
	
	called in calcScore with parameters R[i][0],R[i][1],R[i][2],R[j][0],R[j][1],R[j][2]
		with return variable dist
	called in randCord with parameters rr[i][0],rr[i][1],rr[i][2],0,0,0) and return 
		variable len
	called in calculateCostFunction with parameters temp[j][0],temp[j][1],temp[j][2],R[k][0],R[k][1],R[k][2]
		and return variable rik
	called in calcuAdaptatinScore with parameters R[i][0],R[i][1],R[i][2],temp[j][0],temp[j][1],temp[j][2]
		and return variable disk
*/
double distance(double x1, double y1, double z1, double x2, double y2, double z2)
{
	return (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2);
}




/*
	function to test equality of Probability used in simulated annealing, if P[0~N] are all identical, the vector can choose a random P[?]
	
	no parameters and returns boolean status

	called in growth and used in if statement
	called in adaptation and used in if statement (commented out)
*/
bool isProbSame()
{
	for(int j=0;j<Nt;j++)
		for(int k=j+1;k<Nt;k++)
			if(P[j]!=P[k])
				return false;
	return true;
}




/*
	function to generate random coordinate, used in growth step to generate vectors
	
	no parameters and no return
	
	called in genetic 
	called in main (commented out)
	called in adaptation 
*/
void randCord()
{
	double len, sum = 0.0;
//	srand ( time(NULL) );
	for( int j=0;j<3;j++ ){

		R[0][j] = ((double)rand() - RAND_MAX)*2/RAND_MAX;

	}


	for( int i=0;i<Nt;i++ ){
		for( int j=0;j<3;j++ ){
			rr[i][j] = (double)rand()/RAND_MAX;
		}
		len = distance(rr[i][0],rr[i][1],rr[i][2],0,0,0);
		sum += len;
	}

}



/*
	function to write 3D coordinates info to out pdb file
	
	parameter is pdb file

	no return
	
	called in calcScore with f which is contact info
	called in adaptation with parameter f (commented out)
*/
void outPDB( char* pdbfile )
{
	ofstream pdb;
	pdb.open( pdbfile );
	pdb.precision(3);
	pdb.setf(ios::fixed, ios::floatfield);
	int l,f=3;
	char s[3];
	pdb << "COMPND    NUMBER OF UNITS\t" << n << endl ;
	for( int i=0;i<n;i++ )
	{
		sprintf(s,"%d",i+1);
		l = f-strlen(s);

		pdb << "ATOM" << setw(7) << i+1 << setw(5) << "CA" << setw(4) << i/10+1 << setw(2) << "A" << (i+1) << setw(l+13);
		pdb << R[i][0] << setw(8) << R[i][1] << setw(8) << R[i][2];
		pdb  << setw(6) << "1.00"  << setw(6) << "75.00" << endl;
	}
	for( int i=1;i<n;i++ )
	{
			pdb << "CONECT"<< setw(5) << i << setw(5) << i+1 << endl;
	}
	pdb << "END" << endl;
	pdb.close();
}




/*
	function to calculate the energy function for simulated annealing in growth step

	parameters int i and int j
	
	returns int sum

	called in calculateProbability with parameters i, l and return assigned to Eg
*/
int calculateCostFunction(int i, int j)
{
	int d; int sum = 0;
	double rik;
	for( int k = 1; k <= i-1; k++ )
	{
		d = i-k;
		rik = distance(temp[j][0],temp[j][1],temp[j][2],R[k][0],R[k][1],R[k][2]);		//because temp[] contains rij
		sum += d * aG(S[i][k]) * heaviside(dt-rik);		//d * ag(S[i][k]) * v(dt-r[i][k]);
	}
	return sum;
}




/*
	function to calculate acceptance probability for simulated annealing in growth step

	parameters int i and int j
	
	returns double value of e/sum

	called in growth with parameters i+1, j and return assigned to P[j]
	called in parameters with same parameters and return (commented out)
	called in calcWeight with paramaters i, j and return to double p
*/
double calculateProbability(int i, int j)
{
	int Eg;
	double e, sum = 0.0;
	for( int l = 0; l < Nt; l++ )
	{
		Eg = calculateCostFunction(i, l);	//cout << "Eg=" << Eg << endl;
		e = -Eg/Tg;
		e = exp(e);			//cout << "e=" << e << endl;
		sum += e;			//cout << "sum=" << sum << endl;

	}
	//cout << "sum prob=" << sum << endl;
	Eg = calculateCostFunction(i, j);
	e = -Eg/Tg;
	e = exp(e);				//cout << "e=" << e << endl;
	return e/sum;
}




/*
	function to translate coordinate of the rest of structure accordingly after a unit change its coordinates
	
	parameters int i and l, no return

	called in adaptation with parameters i,j 
*/
void translateCoordinate( int i, int l )
{
//	double t[3];
	for( int k=i+2;k<n;k++ )
	{
		for( int j=0;j<3;j++ )
		{
			R[k][j] = temp[l][j] + (R[k][j]-R[i+1][j]);
		}
	}
}



/*
	function to calculate the score of model for adaptation step

	parameters int I and return int

	called in adaptation with parameters i+1 return value scoreJ
*/
int calcuAdaptatinScore( int I )
{
	double dist;

	int score[Nt]={0};
	int max,maxj;
	max=maxj=0;

	for( int i=0; i<n; i++ )
	{
		if(i==I)		//diagonal
			continue;
		for( int j=0; j<Nt; j++ )
		{
			dist = distance(R[i][0],R[i][1],R[i][2],temp[j][0],temp[j][1],temp[j][2]);
			if(dist<dt && S[i][I])				//reward contact satisfied
					score[j]++;
			else if(dist>dt && !S[i][I])
				if(dist<mymaxd)					//max distance threshold is implemented to make sure the distance is valid
					score[j]++;					//reward non-contact satisfied
			else
				score[j]--;						//penalized unsatisfied contact and non-contact

		}
	}
	for( int j=0; j<Nt; j++ )
	{
		if(score[j]>max)
		{
			max = score[j];
			maxj = j;
		}
	}
	if(max==0)
		return -1;
	else
		return maxj;
}



/*
	function to calculate score and output the coordinates to pdb file if the score of model higher than threshold
	
	parameters: the paramater used here to specify which stage and how many iteration is the pdb file generated
	
	returns double
	
	called multiple times in main with no parameter and no assigning to return(commented out)
	called in genetic with no parameters and no assigning to return
	called in adaptation with no parameters and no return 
		called again adaptation with return value snew
		called again with parameter 1 and no return
*/
void calcScore(int s)
{
	int adapt=0;
	double dist;                  //gets value from distance
	char f[20];                   //used to store info, sent to outpdb
	int count,noncount,contactCount,noncontactCount,non,maxx, tot;
	count=noncount=contactCount=noncontactCount=non=maxx=tot=0;
	int c, pif=0;

	int *pos = new int[n];

	//resetting T[][]
	for(int i=0;i<n;i++)
		for(int j=0;j<i;j++)
			T[i][j] = 0;

	//value DBL_MAX is used, part of library and value is 1.7976931348623157E+308, represents max value of a double
	double mind=DBL_MAX, maxd=0.0, minc=DBL_MAX, maxc=0.0, sum=0.0, nonsum=0.0, uncon=0.0, unnoncon=0.0;
	
	//traverse through one half of matrix S
	for(int i=0;i<n;i++)
	{
		for(int j=0;j<i;j++)
		{
			++tot;
			dist = distance(R[i][0],R[i][1],R[i][2],R[j][0],R[j][1],R[j][2]);
	
			//find max of all the dist and set to maxd
			if(maxd<dist)
				maxd=dist;
			//sum+=dist;
			sum+=dist*S[i][j];

			//now find the min & max dist. of conseq. coords
			if(j==i-1)
			{
//				cout << " [" << i << "][" << j << "]=" << dist <<endl;
				if(minc>dist)
					minc=dist;
				if(maxc<dist)
					maxc=dist;
			}

			//mymaxd is global variable, initially set to 10
			if(dist<=mymaxd)
				maxx++;

			//if contact
			if(S[i][j]){
				contactCount++;
				if(dist<dt){               //distance threshold dt is used and is a global variable, initially set to 5
					count++;
					c=1;
					T[i][j] = 1;
					pif += S[i][j];
				}
				else{
					uncon+=dist;
					T[i][j] = 2;
				}
			}

			//if non contact
			if(!S[i][j]){
				noncontactCount++;
				nonsum+=dist;
				if(dist>dt){                //distance threshold dt is used and is a global variable, initially set to 5
					noncount++;
					c=0;
					T[i][j] = 3;
				}
				else{
					unnoncon+=dist;
					T[i][j] = 4;
				}
			}


		}
	}

	double percent1 = (count/(double)contactCount)*100.0;
	double percent2 = (noncount/(double)noncontactCount)*100.0;
	double avg = sum/totalIF;
	double nonavg = nonsum/noncontactCount;
	double maxdd = (maxx/(double)tot)*100.0;

	double pavg = (pif/(double)totalIF) * 100;

	double unconavg = uncon/(contactCount-count);
	double unnonavg = unnoncon/(noncontactCount-noncount);

	if(maxcon<=percent1 || maxnon<=percent2){
		maxcon=percent1;
		maxnon=percent2;
	}

	if(adapt){

	//write info to log (global variable)
	mylog << "maxd=" << maxd << ",maxdd=" << maxdd << ",maxc=" << maxc << ",avg=" << avg << endl;
	mylog << "actual=" << contactCount << ",real=" << count << "," << percent1 << "%\n";
	mylog << "noncon=" << noncontactCount << ",noncount=" << noncount << "," << percent2 << "%\n";
	mylog << "un con avg=" << unconavg << ",un non-con avg=" << unnonavg << endl;
	mylog << "% of IF satisfied = " << pavg << ", total IF = " << totalIF << ",pif = " << pif <<endl<<endl;
	sprintf(f,"%sF_%.0f_%.0f_%.0f_%.0f.pdb",filename,percent1,percent2,maxdd,pavg);
	if(percent1>=CON+5 || percent2>=NON+5 )
	return;
	}


	if(percent1>CON && percent2>NON )
	{
		sprintf(f,"%sF_%.0f_%.0f_%.0f_%.0f_%i.pdb",filename,percent1,percent2,maxdd,pavg,s);
		outPDB(f);
		mylog << "maxd=" << maxd << ",maxdd=" << maxdd << ",maxc=" << maxc << ",avg=" << avg << endl;
		mylog << "actual=" << contactCount << ",real=" << count << "," << percent1 << "%\n";
		mylog << "noncon=" << noncontactCount << ",noncount=" << noncount << "," << percent2 << "%\n";
		mylog << "un con avg=" << unconavg << ",un non-con avg=" << unnonavg << endl;
		mylog << "% of IF satisfied = " << pavg << ", total IF = " << totalIF << ",pif = " << pif <<endl<<endl;

	}
}



/*
	function to perform adaptation
	
	parameters default parameter int a = -1
	
	no return

	called in calcScore with parameters i and j
	called in main with no parameters
*/
void adaptation(int a=-1)
{
	int i;
	float maxP, minP;
	int maxJ, minJ, randJ, scoreJ;
	int *J;

	for(s=0;s<step;s++)
	{
		maxP=0.0;
		maxJ=0;
		minP=DBL_MAX;
		minJ=0;
		if(a==-1)
			i = rand()%(n-1);			//cout<<i<<endl;	//rand from 0 to n-2
		else
			i = a;

		if(flag)						//reverse the flag for each round
			flag--;
		else
			flag=addsub;

		randCord();

		for(int j=0;j<Nt;j++)
		{
			if(flag)
			{
//				cout <<"addition\n";
				temp[j][0] = R[i][0] + rr[j][0];
				temp[j][1] = R[i][1] + rr[j][1];
				temp[j][2] = R[i][2] + rr[j][2];
			}
			else
			{
//				cout <<"subtraction\n";
				temp[j][0] = R[i][0] - rr[j][0];
				temp[j][1] = R[i][1] - rr[j][1];
				temp[j][2] = R[i][2] - rr[j][2];
			}

		}

		scoreJ = calcuAdaptatinScore(i+1);
		if(scoreJ==-1)
		{
			continue;
		}
		else
		{
			translateCoordinate(i,scoreJ);
			R[i+1][0] = temp[scoreJ][0];
			R[i+1][1] = temp[scoreJ][1];
			R[i+1][2] = temp[scoreJ][2];
		}

		calcScore(s);

	}
	cout<<"adaptation finished\n";
}


//the growth step for chromosome initialization with longer length (>200, chromosome 1,2,3)
void growth()
{
	float maxP;int maxJ, randJ;
	cout<<"growth started\n";
	randCord();
	for(int i=0;i<n-1;i++)
	{
		maxP=0.0;
		for(int j=0;j<Nt;j++)
		{
			temp[j][0] = R[i][0] + rr[j][0];
			temp[j][1] = R[i][1] + rr[j][1];
			temp[j][2] = R[i][2] + rr[j][2];
		}
		for(int j=0;j<Nt;j++)
		{
			P[j] = calculateProbability(i+1, j);	//bcz we r adding (i+1)th monomer

			if( P[j]>maxP )
			{
				maxP = P[j];
				maxJ = j;
			}
		}
		if(isProbSame())         //if all the p are same
		{
			randJ = rand()%Nt;
			R[i+1][0] = temp[randJ][0];
			R[i+1][1] = temp[randJ][1];
			R[i+1][2] = temp[randJ][2];
		}
		else
		{

			R[i+1][0] = temp[maxJ][0];
			R[i+1][1] = temp[maxJ][1];
			R[i+1][2] = temp[maxJ][2];
		}
	}
	cout<<"growth finished\n";
}



/*
	function to initialization for chromosome with shorter length (<200)

	no parameters and no return

	called in main
*/
void initSphereConform()
{
	double x1, x2, x, y, z, x1x2;
	int i = 0;
	while(i!=n)
	{
		x1 = -1+2*((double)rand())/RAND_MAX;
		x2 = -1+2*((double)rand())/RAND_MAX;
		x1x2 = x1*x1 + x2*x2;

		if( (x1x2) >= 1 )
			continue;
		x = 2*x1*sqrt(1-x1x2);
		y = 2*x2*sqrt(1-x1x2);
		z = 1 - 2*(x1x2);
		R[i][0]=x,R[i][1]=y,R[i][2]=z,++i;
	}

}



/*
 * genetic algorithm used to perform crossover and mutation
 *
*/
void genetic()
{
	int cp = rand()%n;
	double temp[n][3];
	cout << cp <<endl;

	//print values of A and B
	for(int i=0;i<n;i++) {
		cout << A[i][0] << " " << A[i][1] << " " << A[i][2] << "\t" << B[i][0] << " " << B[i][1] << " " << B[i][2] << endl;
	}
	//assign A to temp, B to A, and temp to B, mutation step
	for( int i=0;i<n;i++ ){
		temp[i][0] = A[i][0];
		temp[i][1] = A[i][1];
		temp[i][2] = A[i][2];
	}

	for( int i=cp;i<n;i++ ){
		A[i][0] = B[i][0];
		A[i][1] = B[i][1];
		A[i][2] = B[i][2];
	}
	for( int i=cp;i<n;i++ ){
		B[i][0] = temp[i][0];
		B[i][1] = temp[i][1];
		B[i][2] = temp[i][2];
	}
	//print A and B again
	for(int i=0;i<n;i++) {
		cout << A[i][0] << " " << A[i][1] << " " << A[i][2] << "\t" << B[i][0] << " " << B[i][1] << " " << B[i][2] << endl;
	}
	//assign A to R
	for( int i=0;i<n;i++ ){
		R[i][0] = A[i][0];
		R[i][1] = A[i][1];
		R[i][2] = A[i][2];
	}
	calcScore(-3);
	//assign B to R
	for( int i=0;i<n;i++ ){
		R[i][0] = B[i][0];
		R[i][1] = B[i][1];
		R[i][2] = B[i][2];
	}
	calcScore(-4);

	//mutation: randomly selection element in A and B and change vale to that of rr
	randCord();
	int m = rand()%n;
	cout << m <<endl;

	A[m][0] += rr[0][0];
	A[m][1] += rr[0][1];
	A[m][2] += rr[0][2];

	m = rand()%n;
	cout << m <<endl;
	B[m][0] += rr[9][0];
	B[m][1] += rr[9][1];
	B[m][2] += rr[9][2];

	//print A and B
	for(int i=0;i<n;i++) {
		cout << A[i][0] << " " << A[i][1] << " " << A[i][2] << "\t" << B[i][0] << " " << B[i][1] << " " << B[i][2] << endl;
	}

}




//function used to input the 3D strcture in the pdb file, used for rubustness test
void inPDB( char* pdbfile )
{
	ifstream myfile;
	string line;

	char *result=NULL;
	char *ln;
	int i = 0;

	myfile.open( pdbfile );
	while( !myfile.eof() )
	{
		getline(myfile, line);
		ln = (char *) line.c_str();
		result = strtok(ln," ");

		while( result != NULL ) {
			if(!strcmp(result,"ATOM"))
			{
				result = strtok( NULL, " " );
				result = strtok( NULL, " " );
				result = strtok( NULL, " " );
				result = strtok( NULL, " " );
				R[i][0] = atof(strtok( NULL, " " ));
				R[i][1] = atof(strtok( NULL, " " ));
				R[i][2] = atof(strtok( NULL, " " ));

				cout << R[i][0] << " " << R[i][1] << " " << R[i][2] << endl;
				++i;
			}
			result = strtok( NULL, " " );
		}
	}
	myfile.close();
}




/*
	function that specifies which regions to omit 

	because they are centromere regions, this 

	depends on the chromosome number
	
	parameters include m, value assigned in main 
	based on i, which is chromosome number taken 
	from input file. m is assigned to value mid.
	
	no return value

	called in main
*/
void omitRegion(int m)
{
	int** SS; int mid=m;              //creates a matrix SS size n-oz	
	SS = new int* [n-oz];
	for(int i=0;i<n-oz;i++)
		SS[i] = new int[n-oz];

// initialize SS based on omitting top, left region or mid region
// this is for top, left omit
	if(!m){
	for( int i=oz;i<n;i++ )
		for( int j=oz;j<n;j++ )
			SS[i-oz][j-oz]=S[i][j];
	}
// this is for top, left omit

	else{

 // this is for mid region omit
	for( int i=0,ii=0;i<n;i++ ){
		for( int j=0,jj=0;j<n;j++ ){
			if(j==mid){
				j=j+oz-1;
				continue;
			}
			SS[ii][jj++]=S[i][j];

		}

		if(i==mid){
			i=i+oz-1;
			continue;
		}
		++ii;
	}
// this is for mid region omit
	}

	//create two more matrices, S with n-oz rows and n-oz columns and
	//							R with n-oz rows and 3 columns
	S = new int* [n-oz];
	R = new double* [n-oz];
	for(int i=0;i<n-oz;i++){
		S[i] = new int[n-oz];
		R[i] = new double[3];
	}

	n=n-oz;					         //changing n to omit regions
	for( int i=0;i<n;i++ )           //copy values from SS to S
		for( int j=0;j<n;j++ )
			S[i][j]=SS[i][j];

	//delete all of matrix SS
	for(int i=0;i<n;i++)
		delete []SS[i];
	delete []SS;

}





int main( int argc, char *argv[] )
{
	// type ofstream
	// precision used to set decimal places of double

	if( argc != 2 )
	{
		cout<<"usage: "<<argv[0]<<" <interaction_file>"<<endl;
		return 1;
	}


/*
	if( argc != 3 )				//used for rebustness test, the third paramater can be a pdb file which you can used it as initial structure to generate models
	{
		cout<<"usage: "<<argv[0]<<" <interaction_file> <pdb_file>"<<endl;
		return 1;
	}
*/
	mylog.open( "score.log" );
	mylog.precision(3);


	cout << "hello world"<< endl;
	srand ( time(NULL) );

	// value later assigned, depends on value of i
	// sent to omitRegion
	int m;

	// input file
	// type char*
	filename = argv[1];
	cout<<filename<<endl;

	int i = atoi(filename);               // assigns i to integer value of filename (input file)

	mylog << "chromosome=" <<i<<endl;

	// read the input to figure out number of elements for matrices and
	// assign data from input to said matrices
	readInput();

	// based on value of i, assign oz and m accordingly and 
	// send the right value to omitRegion
	if(i==1){
		oz=19;
		m=122;
		omitRegion(m);
	}
	else if(i==2||i==3){
		oz=i;
		m=94-i;
		omitRegion(m);
	}
	else if(i==4){
		oz=2;
		m=49;
		omitRegion(m);
	}
	else if(i==5){
		oz=2;
		m=47;
		omitRegion(m);
	}
	else if(i==6||i==7){
		oz=2;
		m=59;
		omitRegion(m);
	}
	else if(i==8){
		oz=2;
		m=44;
		omitRegion(m);
	}
	else if(i==9){
		oz=17;
		m=48;
		omitRegion(m);
	}
	else if(i==10){
		oz=1;
		m=40;
		omitRegion(m);
	}
	else if(i==11){
		oz=2;
		m=52;
		omitRegion(m);
		}
	else if(i==12){
		oz=1;
		m=36;
		omitRegion(m);
	}
	else if(i==13||i==14){
		oz=i+4;
		m=0;
		omitRegion(m);
	}
	else if(i==15){
		oz=18;
		m=0;
		omitRegion(m);
	}
	else if(i==16){
		oz=8;
		m=36;
		omitRegion(m);
	}
	else if(i==19){
		oz=7;
		m=25;
		omitRegion(m);
	}
	else if(i==20){
		oz=1;
		m=27;
		omitRegion(m);
	}
	else if(i==21){
		oz=18;
		m=0;
		omitRegion(m);
		oz=2;m=2;
		omitRegion(m);
	}
	else if(i==22){
		oz=14;
		m=0;
		omitRegion(m);
	}

	printContactMatrix();

	if(i==1||i==2||i==3)	//base on choromsome length > or < 200? to choose the intialization method, chromosome 1,2,3's length is longer than 200
		growth();
	else
		initSphereConform();

/*
	inPDB("14F_65_72_85_80_24075.pdb");			//input you strcture for genetic algorithm for genetic parents
	printR(0);
	calcScore(-5);
	inPDB("14F_65_72_85_80_24077.pdb");
	printR(1);
	calcScore(-6);
	genetic();
*/

/*
	inPDB(argv[2]);		//used for rebustness test, the third paramater can be a pdb file which you can used it as initial structure to generate models,
						//in this case, the growth or Sphere is not needed since already have the pdb as initial structure,
*/
	calcScore(-1);

	adaptation();		//adaptation step, which can be command out when using genetic algorithm
	
	calcScore(-2);

	cout<<"maxcon="<<maxcon<<",maxnon="<<maxnon<<endl;

	memdel();                     //delete matrices and their data
	mylog.close();                //close log

	cout << "end..." << endl;
	return 0;
}
