#include <iostream>
#include <cstdio>
#include <fstream>
#include <ctime>
#include <random>

/******
        16/12/19

        Parameter extraction for solar cell operation using Genetic Algorithm using C++.
        5 parameters in the single diode equation: Il, I0, eta, Rs, Rsh.
        by Prachi Bisht

        (see pdf for theory)
*******/

std::mt19937 gen{static_cast<long unsigned int>(time(0))};
std::uniform_real_distribution<double> dist_u(0,1);
std::uniform_int_distribution<int> dist_80(0,79);
std::uniform_int_distribution<int> dist_1000(0,999);

const double q = 1.6e-19;               //fundamental charge
const double kb = 1.38e-23;             //Boltzmann constant
const double T = 300.0;                 //room temperature

const int n = 1000;                     //the size of gene pool

int top_genome_no = 333;                //fraction of gene pool used for reproducing next generation of genes
int elite_no = 5;                       //number of genes going straight to the next gen
int generation_no = 500;                //number of iterations. increase this for more accuracy

/*upper and lower limits for the five parameters. Intervals should contain the solution.
  Narrower the interval, higher the accuracy*/
double ilmin = 0.400;
double ilmax = 0.600;
double i0min = 1e-14;
double i0max = 1e-9;
double nmin = 1;
double nmax = 3;
double rsmin = 0.0;
double rsmax = 0.1;
double rshmin = 10000.0;
double rshmax = 1000000.0;

/*here 80 = 16*5 represents the array of bits representing one gene, 16 bits for each of the 5 traits (Il, I0...)*/
int genome[n][80], child[n][80];
double il_m[n], i0_m[n], n_m[n], rs_m[n], rsh_m[n];
double err_m[n];

/*initialize the genome randomly*/
void initialize (int (*genome)[80]);

/*find the values of the parameters encoded by the gene in decimal*/
void map_values( int (*genome)[80], double *il_m, double *i0_m, double *n_m, double *rs_m, double *rsh_m);

/*difference between measured current and current calculated using the parameters encoded by a gene*/
void calc_error(double *err_m, double *voltage, double *current, double *il_m, double *i0_m, double *n_m, double *rs_m, double *rsh_m );

/*sort the genes according to their fitness*/
void sort(double *err_m, int (*genome)[80]);

/*transfer top 5 fittest genes to the next iteration and populate the rest by randomly crossing between
top 30 % of the current iteration*/
void populate(int (*genome)[80], int (*child)[80]);

/*randomly mutate some bits once in a while; allow for stochastic sampling*/
void mutate(int (*child)[80]);


int main()
{
	/*store VIdata. More granular the data, better the accuracy.*/

	int m = 0;

    /*data from experiment/TCAD simulation*/
    std::ifstream fin;
    fin.open("gaas.txt",std::ios::in);

    std::string line;
    for(m = 0; std::getline(fin,line); m++);
    fin.clear();
    fin.seekg(0, std::ios::beg);

    double current[m], voltage[m];

	long double t = 0.0;
	int j = 0; int k = 0;


    while(j<2*m)
	{
		fin>>t;
		if (j % 2 == 0)
			voltage[k] = t;
		else
		{
			current[k] = t;
			++k;
		}
		j++;

	}
	fin.close();

	/*initialize the genome set with random bit values for the variables*/
	initialize(genome);

	for(int gen = 0; gen<generation_no; ++gen)
	{
		/*finding values for all five variables*/
	   map_values(genome, il_m, i0_m, n_m, rs_m, rsh_m);


		/*calculating error*/
	   calc_error(err_m, voltage, current, il_m, i0_m, n_m, rs_m, rsh_m);


		/*sorting error values in increasing order. best genomes are on top*/
	   sort(err_m, genome);


		/*retain top 5 of the good genome, the elite ones*/
		int elite[elite_no][80];
		for (int i = 0; i<elite_no; ++i)
		   	for(int j = 0; j<80; ++j)
		   		elite[i][j] = genome[i][j];



	   /*cross the rest*/
	   populate(genome, child);


	   /*mutate once in a while*/
	   mutate(child);
	   mutate(child);


	   /*new generation*/
	   for (int i = 0; i<elite_no; ++i)
		   	for(int j = 0; j<80; ++j)
			   genome[i][j] = elite[i][j];

		for (int i = elite_no; i<n; ++i)
		   	for(int j = 0; j<80; ++j)
		   	   genome[i][j] = child[i][j];


		/*rinse. lather. repeat*/

    }
		double fil = 0.0;
		for( int j = 0; j<16; ++j)
			fil = 2.0*fil + genome[0][j];


		double fi0 = 0.0;
		for( int j = 16; j<32; ++j)
			fi0 = 2.0*fi0 + genome[0][j];


		double fn = 0.0;
		for( int j = 32; j<48; ++j)
			fn = 2.0*fn + genome[0][j];


		double frs = 0.0;
		for( int j = 48; j<64; ++j)
			frs = 2.0*frs + genome[0][j];


		double frsh = 0.0;
		for( int j = 64; j<80; ++j)
			frsh = 2.0*frsh + genome[0][j];


		fil = ilmin + (ilmax - ilmin)*fil/pow(2.0,16);
		fi0 = i0min + (i0max - i0min)*fi0/pow(2.0,16);
		fn = nmin + (nmax - nmin)*fn/pow(2.0,16);
		frs = rsmin + (rsmax - rsmin)*frs/pow(2.0,16);
		frsh = rshmin + (rshmax - rshmin)*frsh/pow(2.0,16);


		double errt = 0.0;

		for (int j = 0; j<m; ++j)
		{

			double vd =(voltage[j] + current[j]*frs);
			double term1 = fil;
			double term2 = fi0*(exp(vd*q/(fn*kb*T))-1);
			double term3 = vd/(1.0*frsh);

			double error = pow(((current[j] - term1 + term2 + term3)/(1.0*current[j])), 2);
			errt += error;
		}

    std::cout<<"After "<<generation_no<<" iterations:"<<"\n";
    std::cout<<"IL = "<<" "<<fil<<" A\n"
            <<"I0 = "<<" "<<fi0<<" A\n"
            <<"eta = "<<" "<<fn<<"\n"
            <<"Rs = "<<" "<<frs<<" ohm\n"
            <<"Rsh = "<<" "<<frsh<<" ohm \n";

	//std::cout<<fil<<" "<<fi0<<" "<<fn<<" "<<frs<<" "<<frsh;
	std::cout<<"Error:"<<" "<<sqrt(errt/(1.0*m))<<"%\n";

	return 0;
}

/*initialize the genome randomly*/
void initialize (int genome[n][80])
{
	for( int i = 0; i<n; ++i)
	{
		for (int j = 0; j<80; ++j)
		{
			double bitvalue = dist_u(gen) + 0.5;
			genome[i][j] = (int) bitvalue;
		}
	}
}


void map_values( int genome[n][80], double *il_m, double *i0_m, double *n_m, double *rs_m, double *rsh_m)
{
	for ( int i = 0; i<n; ++i)
	{
		double il = 0.0;
		for( int j = 0; j<16; ++j)
			il = 2.0*il + genome[i][j];


		double i0 = 0.0;
		for( int j = 16; j<32; ++j)
			i0 = 2.0*i0 + genome[i][j];


		double n = 0.0;
		for( int j = 32; j<48; ++j)
			n = 2.0*n + genome[i][j];


		double rs = 0.0;
		for( int j = 48; j<64; ++j)
			rs = 2.0*rs + genome[i][j];


		double rsh = 0.0;
		for( int j = 64; j<80; ++j)
			rsh = 2.0*rsh + genome[i][j];


		/*mapping bitvalue to actual values b/w the specified range*/


		il = ilmin + (ilmax - ilmin)*il/pow(2.0,16);
		i0 = i0min + (i0max - i0min)*i0/pow(2.0,16);
		n = nmin + (nmax - nmin)*n/pow(2.0,16);
		rs = rsmin + (rsmax - rsmin)*rs/pow(2.0,16);
		rsh = rshmin + (rshmax - rshmin)*rsh/pow(2.0,16);

		/*store the actual values*/

		il_m[i] = il;
		i0_m[i] = i0;
		n_m[i] = n;
		rs_m[i] = rs;
		rsh_m[i] = rsh;


	}
}


void calc_error(double *err_m, double *voltage, double *current, double *il_m, double *i0_m, double *n_m, double *rs_m, double *rsh_m )
{
   	int m = sizeof(voltage)/sizeof(voltage[0]);

	for (int i = 0; i<n; ++i)
	{
		double errt = 0.0;

		for (int j = 0; j<m; ++j)
		{

			double vd = (voltage[j] + current[j]*rs_m[i]);
			double term1 = il_m[i];
			double term2 = i0_m[i]*(exp(vd*q/(n_m[i]*kb*T))-1);
			double term3 = vd/(1.0*rsh_m[i]);
            /*square error*/
			double error = pow(((current[j] - term1 + term2 + term3)/(1.0*current[j])), 2);
			errt += error;
		}

		err_m[i] = sqrt(errt);
	}
}


void sort(double *err_m, int genome[n][80])
{
	int counter = 0;

    do{
    	counter = 0;

	    for (int i = 0; i<n-1; ++i)
	    {
	    	if(err_m[i]>err_m[i+1])
	    	{
	    		double a = err_m[i];
	    		err_m[i] = err_m[i+1];
	    		err_m[i+1] = a;

	    		for(int bit = 0; bit < 80; ++bit)
	    		{
	    			int b = genome[i][bit];
					genome[i][bit] = genome[i+1][bit];
					genome[i+1][bit] = b;
				}

	    		counter++;
			}
		}
	}while(counter!=0);


}


void populate(int genome[n][80], int child[n][80])
{
	for (int i = 0; i<n-1; i = i+2)
	{
		int father_no = 0;
		int mother_no = 0;

		do{
			father_no = dist_1000(gen);
		}while(father_no >= top_genome_no);

		do{
			mother_no = dist_1000(gen);
		}while(mother_no >= top_genome_no);


		int parent_frac = dist_80(gen);

	    for( int j = 0; j < parent_frac; ++j)
		{
			child[i][j] = genome[father_no][j];
			child[i+1][j] = genome[mother_no][j];
		}

		for( int j = parent_frac; j < 80; ++j)
		{
			child[i][j] = genome[mother_no][j];
			child[i+1][j] = genome[father_no][j];
		}

    }

}


void mutate(int child[n][80])
{
	for(int i = 0; i<n; ++i)
	{
		for(int j = 0; j<80; ++j)
		{
			double eta = dist_u(gen);
			if(eta<0.05)
				if(child[i][j]==0)
					child[i][j]=1;
				else
					child[i][j]=0;
		}
	}
}

