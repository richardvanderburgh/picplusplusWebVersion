#include <gtest/gtest.h>
#include <init.hpp>

// Define a test case
//TEST(GTestExample, AdditionTest) {
//    // Test assertion
//    
//    EXPECT_EQ(2 + 2, 4);
//}

//TEST(FFTTest, fftTest) {
//
//	const int ng = 128;
//	std::vector<double> rhoRE(ng, 0.0);
//	std::vector<double> rhoCO(ng, 0.0);
//	std::vector<double> rhokRE(ng, 0.0);
//	std::vector<double> rhokCO(ng, 0.0);
//
//	std::vector<double> phiRE(ng, 0.0);
//	std::vector<double> phiCO(ng, 0.0);
//
//	std::vector<double> space(ng, 0.0);
//
//	space = linspace(0, 2*M_PI, ng);
//
//	double val = sin(55);
//
//	complex complexRho[ng];
//	complex complexRhok[ng];
//	complex complexPhik[ng];
//	complex complexPhi[ng];
//
//	for (int i = 0; i < ng; i++) {
//		complexRho[i] = sin(space[i]);
//		rhoRE[i] = complexRho[i].re();
//		rhoCO[i] = complexRho[i].im();
//
//		std::cout << "rhoRE[" << i << "] = " << rhoRE[i] << std::endl;
//		std::cout << "rhoCO[" << i << "] = " << rhoCO[i] << std::endl;
//	}
//
//	matplot::scatter(space, rhoRE, 1);
//	matplot::scatter(space, rhoCO, 1);
//
//	CFFT::Forward(complexRho, complexRhok, ng);
//
//	for (int i = 0; i < ng; i++) {
//		//complexRho[i] = sin(space[i]);
//		rhokRE[i] = complexRhok[i].re();
//		rhokCO[i] = complexRhok[i].im();
//
//		std::cout << "rhokRE[" << i << "] = " << rhokRE[i] << std::endl;
//		std::cout << "rhokCO[" << i << "] = " << rhokCO[i] << std::endl;
//	}
//
//	matplot::scatter(space, rhokRE, 1);
//	matplot::scatter(space, rhokCO, 1);
//
//	CFFT::Inverse(complexRhok, complexPhi, ng);
//
//	for (int i = 0; i < ng; i++) {
//		//complexRho[i] = sin(space[i]);
//		phiRE[i] = complexPhi[i].re();
//		phiCO[i] = complexPhi[i].im();
//
//		std::cout << "phiRE[" << i << "] = " << phiRE[i] << std::endl;
//		std::cout << "phiCO[" << i << "] = " << phiCO[i] << std::endl;
//	}
//
//	matplot::scatter(space, phiRE, 1);
//	matplot::scatter(space, phiCO, 1);
//}
//
TEST(PICTest, EFrame0Test)
{
	Init init;

	int N = 5;
	int nt = 3;
	double dt = 0.1;

	int mode = 1;
	int V0 = 1;
	int numSpecies = 2;
	double amplitude = 0.001;
	double VT1 = 0;

	bool success = init.initialize(N, nt, dt, mode, V0, numSpecies, amplitude, VT1);

	int ng = 32;

	std::vector<double> expectedE{
	   0.00004269726342988,
	   0.40061715226629907, 
	   0.69129615923052634, 
	   0.15270125999750908, 
	 - 0.58858351024757993, 
	 - 0.53469273983236021, 
	 - 0.15700368113596019, 
	   0.23764672981919918, 
	   0.5872089012479681 ,
	   0.44276831815030276, 
	 - 0.29787742601334988, 
	 - 0.63964937826384105, 
	 - 0.31838214086949496, 
	   0.07641105814515342,
	   0.48233137539586241, 
	   0.73622824655844488, 
	 - 0.00577633017558041,
	 - 0.74269508107571602, 
	 - 0.4829274975673129 ,
	 - 0.07637179844676961,
	   0.31841567826749095, 
	   0.63923533770152829, 
	   0.29429815402774046, 
	 - 0.44632825059627845, 
	 - 0.58765143316935486, 
	 - 0.23755375359532441, 
	   0.15685674017510517, 
	   0.5358055499237665 ,
	   0.59793227836901386, 
	 - 0.14335563303744897, 
	 - 0.69017996479850396, 
	 - 0.40076701771446444, 
	   0.00004269726342988
	};

	bool fieldsMatch = true;
	double precision = 10E-12;

	double actualE;
	for (int i = 0; i < ng; i++) {

		actualE = init.mPicData.frames[0].electricField[i];

		if ((actualE - expectedE[i]) > precision) {
			fieldsMatch = false;
		}

	}


	EXPECT_TRUE(fieldsMatch);
}

TEST(PICTest, EFrame1Test)
{
	Init init;

	int N = 5;
	int nt = 3;
	double dt = 0.1;

	int mode = 1;
	int V0 = 1;
	int numSpecies = 2;
	double amplitude = 0.001;
	double VT1 = 0;

	bool success = init.initialize(N, nt, dt, mode, V0, numSpecies, amplitude, VT1);

	int ng = 32;

	std::vector<double> expectedE{
	 -0.0000439013224624,
	  0.3788884122956504,
	  0.5380417214081206,
	  0.1530739948387394,
	-0.4363518996148261,
	-0.5113511230499872,
	-0.1608469453677030,
	  0.2401641546222603,
	  0.5782959113535076,
	  0.3895727295358895,
	-0.2986005920314539,
	-0.5846991593056296,
	-0.3123233273425899,
	  0.0807179914645374,
	  0.4447460941407272,
	  0.4840209388191416,
	-0.0057212846921208,
	-0.4904166609622277,
	-0.4453481847763712,
	-0.0806466975043333,
	  0.3121236745094740,
	  0.5827285304205471,
	  0.2950614898913002,
	-0.3915460838758390,
	-0.5784859556413913,
	-0.2401348554053200,
	  0.1607979225032119,
	  0.5118367972706172,
	  0.4415258488968978,
	-0.1438156987461032,
	-0.5328705719143205,
	-0.3783932704179430,
	  0.0000439013224624
	};

	bool fieldsMatch = true;
	double precision = 10E-12;

	double actualE;
	for (int i = 0; i < ng; i++) {

		actualE = init.mPicData.frames[1].electricField[i];

		if ((actualE - expectedE[i]) > precision) {
			fieldsMatch = false;
		}

	}


	EXPECT_TRUE(fieldsMatch);
}

TEST(PICTest, E_N_5_nt_15_Frame0Test)
{
	Init init;

	int N = 5;
	int nt = 3;
	double dt = 0.1;

	int mode = 1;
	int V0 = 1;
	int numSpecies = 2;
	double amplitude = 0.001;
	double VT1 = 0;

	bool success = init.initialize(N, nt, dt, mode, V0, numSpecies, amplitude, VT1);

	int ng = 32;

	std::vector<double> expectedE{
-0.00007128172188 ,
0.34363433331871  ,
0.25798036035269  ,
0.03511076599531  ,
- 0.17492056009884,
- 0.35437214364456,
- 0.14620560851485,
0.21259948643449  ,
0.31998452983465  ,
0.12698150956722  ,
- 0.08068314890847,
- 0.28786437014873,
- 0.27907308194332,
0.07963121207181  ,
0.38881510288119  ,
0.22078685815613  ,
- 0.00065933311229,
- 0.22625195450655,
- 0.39439842648401,
- 0.08025083508866,
0.27894035077967  ,
0.28618504218963  ,
0.07881315865435  ,
- 0.12894715121960,
- 0.32168365721081,
- 0.21288820391492,
0.14679100015775  ,
0.35883716565894  ,
0.17998365440103  ,
- 0.03012850899719,
- 0.25361305906686,
- 0.34306320587203,
- 0.00007128172188
	};

	bool fieldsMatch = true;
	double precision = 10E-12;

	double actualE;
	for (int i = 0; i < ng; i++) {

		actualE = init.mPicData.frames[2].electricField[i];

		if ((actualE - expectedE[i]) > precision) {
			fieldsMatch = false;
		}

	}


	EXPECT_TRUE(fieldsMatch);
}
