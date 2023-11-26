#include <gtest/gtest.h>

#include <nlohmann/json.hpp>

#include <DataStructs.h>
#include <PICPlusPlus.h>

namespace {
	using DBL_VEC = std::vector<double>;
	using JSON_VEC = std::vector<nlohmann::json>;

	nlohmann::json BuildFrame0Json()
	{
		return nlohmann::json{{"frameNumber", 0},
							  {"electricField", DBL_VEC{4.269726342988564e-05, 0.4006171522662996, 0.6912961592305251, 0.15270125999750908, -0.5885835102475788, -0.5346927398323591, -0.15700368113595964, 0.23764672981919863, 0.5872089012479687, 0.4427683181503034, -0.297877426013351, -0.639649378263841, -0.3183821408694944, 0.07641105814515342, 0.48233137539586185, 0.7362282465584438, -0.005776330175579284, -0.7426950810757155, -0.482927497567314, -0.07637179844676961, 0.31841567826749095, 0.6392353377015277, 0.29429815402774046, -0.446328250596279, -0.5876514331693549, -0.2375537535953244, 0.15685674017510404, 0.535805549923767, 0.5979322783690139, -0.14335563303744897, -0.6901799647985029, -0.400767017714465, 4.269726342988564e-05}},
							  {"particles", JSON_VEC{
												{{"id", 0}, {"position", 3.204120289718404}, {"species", 0}, {"velocity", 0.5093312139007288}},
												{{"id", 1}, {"position", 9.598426189371066}, {"species", 0}, {"velocity", 0.5092842698099325}},
												{{"id", 2}, {"position", 15.994907041821058}, {"species", 0}, {"velocity", 0.5092449562549437}},
												{{"id", 3}, {"position", 22.398426189371065}, {"species", 0}, {"velocity", 0.5092757821427519}},
												{{"id", 4}, {"position", 28.804120289718405}, {"species", 0}, {"velocity", 0.5093428673619713}},
												{{"id", 5}, {"position", 3.204120289718404}, {"species", 1}, {"velocity", -0.5092604218874025}},
												{{"id", 6}, {"position", 9.598426189371066}, {"species", 1}, {"velocity", -0.5093073659781988}},
												{{"id", 7}, {"position", 15.994907041821058}, {"species", 1}, {"velocity", -0.5093466795331876}},
												{{"id", 8}, {"position", 22.398426189371065}, {"species", 1}, {"velocity", -0.5093158536453793}},
												{{"id", 9}, {"position", 28.804120289718405}, {"species", 1}, {"velocity", -0.50924876842616}}}}};
	}

	nlohmann::json BuildFrame1Json()
	{
		return nlohmann::json{{"frameNumber", 1},
							  {"electricField", DBL_VEC{-4.3901322461342085e-05, 0.3788884122956493, 0.538041721408119, 0.15307399483873943, -0.4363518996148244, -0.5113511230499866, -0.16084694536770308, 0.24016415462226148, 0.5782959113535082, 0.389572729535889, -0.29860059203145506, -0.5846991593056291, -0.3123233273425894, 0.08071799146453694, 0.44474609414072663, 0.4840209388191416, -0.0057212846921197635, -0.49041666096222714, -0.4453481847763724, -0.08064669750433388, 0.3121236745094752, 0.5827285304205461, 0.295061489891298, -0.3915460838758385, -0.578485955641389, -0.24013485540531948, 0.16079792250320965, 0.5118367972706173, 0.44152584889690005, -0.14381569874610328, -0.5328705719143222, -0.378393270417943, -4.3901322461342085e-05}},
							  {"particles", JSON_VEC{
												{{"id", 0}, {"position", 3.7133807116058066}, {"species", 0}, {"velocity", 0.5092604218874025}},
												{{"id", 1}, {"position", 10.107733555349265}, {"species", 0}, {"velocity", 0.5093073659781989}},
												{{"id", 2}, {"position", 16.504253721354246}, {"species", 0}, {"velocity", 0.5093466795331877}},
												{{"id", 3}, {"position", 22.907742043016444}, {"species", 0}, {"velocity", 0.5093158536453791}},
												{{"id", 4}, {"position", 29.313369058144566}, {"species", 0}, {"velocity", 0.50924876842616}},
												{{"id", 5}, {"position", 2.6947890758176754}, {"species", 1}, {"velocity", -0.5093312139007288}},
												{{"id", 6}, {"position", 9.089141919561134}, {"species", 1}, {"velocity", -0.5092842698099324}},
												{{"id", 7}, {"position", 15.485662085566114}, {"species", 1}, {"velocity", -0.5092449562549436}},
												{{"id", 8}, {"position", 21.889150407228314}, {"species", 1}, {"velocity", -0.5092757821427522}},
												{{"id", 9}, {"position", 28.294777422356432}, {"species", 1}, {"velocity", -0.5093428673619713}}}}};
	}

	nlohmann::json BuildFrame2Json()
	{
		return nlohmann::json{{"frameNumber", 2},
							  {"electricField", DBL_VEC{-7.128172187758459e-05, 0.3436343333187128, 0.25798036035268507, 0.035110765995305174, -0.17492056009884052, -0.35437214364456043, -0.14620560851484807, 0.21259948643448578, 0.31998452983464865, 0.12698150956722631, -0.08068314890847027, -0.2878643701487316, -0.2790730819433162, 0.07963121207180539, 0.38881510288118526, 0.2207868581561343, -0.0006593331122882585, -0.22625195450654836, -0.39439842648401086, -0.08025083508865681, 0.2789403507796752, 0.2861850421896325, 0.0788131586543521, -0.12894715121960149, -0.3216836572108103, -0.21288820391491628, 0.14679100015775073, 0.35883716565893686, 0.17998365440102762, -0.030128508997190877, -0.25361305906686255, -0.34306320587203326, -7.128172187758459e-05}},
							  {"particles", JSON_VEC{
												{{"id", 0}, {"position", 4.236260267412681}, {"species", 0}, {"velocity", 0.5228795558068746}},
												{{"id", 1}, {"position", 10.633818297349798}, {"species", 0}, {"velocity", 0.5260847420005327}},
												{{"id", 2}, {"position", 17.026339454455755}, {"species", 0}, {"velocity", 0.5220857331015083}},
												{{"id", 3}, {"position", 23.43377303999503}, {"species", 0}, {"velocity", 0.5260309969785861}},
												{{"id", 4}, {"position", 29.83615152044991}, {"species", 0}, {"velocity", 0.5227824623053451}},
												{{"id", 5}, {"position", 2.1716778271241517}, {"species", 1}, {"velocity", -0.5231112486935235}},
												{{"id", 6}, {"position", 8.563141153382492}, {"species", 1}, {"velocity", -0.5260007661786419}},
												{{"id", 7}, {"position", 14.963879707115325}, {"species", 1}, {"velocity", -0.5217823784507886}},
												{{"id", 8}, {"position", 21.363223235794123}, {"species", 1}, {"velocity", -0.5259271714341911}},
												{{"id", 9}, {"position", 27.77173549692073}, {"species", 1}, {"velocity", -0.5230419254357017}}}}};
	}

	nlohmann::json BuildFrame3Json()
	{
		return nlohmann::json{{"frameNumber", 3},
							  {"electricField", DBL_VEC{0.00048611652793800227, 0.12486578690253407, -0.048241322416612546, -0.05152435356781079, 0.12017387971269856, -0.05254602946667631, -0.10615708008743974, 0.17274652702089568, 0.014771533002604689, -0.1553839279452945, 0.10359012009747204, 0.019130542133291885, -0.15427685478175282, 0.05221627969517066, 0.08703164021148174, -0.08589210980149622, -0.0006210014556194892, 0.08069261111846183, -0.09230027776019586, -0.05301134986936989, 0.15278304188257313, -0.020686110679284365, -0.10396259300278209, 0.15369654852543807, -0.01645534177997643, -0.17300322661546264, 0.10668724656003992, 0.05689858912338888, -0.11588310866014291, 0.052524745393571354, 0.05237400194971378, -0.12072452196735774, 0.00048611652793800227}},
							  {"particles", JSON_VEC{
												{{"id", 0}, {"position", 4.7702077298861525}, {"species", 0}, {"velocity", 0.5339474624734711}},
												{{"id", 1}, {"position", 11.170700029090199}, {"species", 0}, {"velocity", 0.5368817317404001}},
												{{"id", 2}, {"position", 17.560173666308433}, {"species", 0}, {"velocity", 0.5338342118526795}},
												{{"id", 3}, {"position", 23.970629173127335}, {"species", 0}, {"velocity", 0.5368561331323038}},
												{{"id", 4}, {"position", 30.369985470545668}, {"species", 0}, {"velocity", 0.5338339500957591}},
												{{"id", 5}, {"position", 1.6373764023950277}, {"species", 1}, {"velocity", -0.534301424729124}},
												{{"id", 6}, {"position", 8.026379140425227}, {"species", 1}, {"velocity", -0.5367620129572658}},
												{{"id", 7}, {"position", 14.430543642998291}, {"species", 1}, {"velocity", -0.5333360641170333}},
												{{"id", 8}, {"position", 20.826556912400264}, {"species", 1}, {"velocity", -0.5366663233938587}},
												{{"id", 9}, {"position", 27.2374478328234}, {"species", 1}, {"velocity", -0.5342876640973317}}}}};
	}

	nlohmann::json BuildExpectedJson()
	{
		return nlohmann::json{
			{"ese", DBL_VEC{0.6490120680288294, 0.4801652664534216, 0.1777764530230646, 0.03045331428993565}},
			{"ke", std::vector<DBL_VEC>{DBL_VEC{0.8148733127700093, 0.8148733127700092, 0.8625254009474299, 0.899446802247833}, DBL_VEC{0.8148733127700093, 0.8148733127700093, 0.8625250675472964, 0.899446068825189}}},
			{"phaseFrames", JSON_VEC{BuildFrame0Json(), BuildFrame1Json(), BuildFrame2Json(), BuildFrame3Json()}}};
	}
}


TEST(PICTest, EFrame0Test)
{
	double spatialLength = 6.28318530717958;
	int numParticles = 5;
	int numTimeSteps = 3;
	double timeStepSize = 0.1;
	int numGrid = 32;
	int spatialPerturbationMode = 1;
	int driftVelocity = 1;
	int numSpecies = 2;
	double spatialPerturbationAmplitude = 0.001;
	double thermalVelocity = 0;
	double plasmaFrequency = 1.0;
	double chargeMassRatio = -1.0;

	std::vector<DATA_STRUCTS::SpeciesData> allSpeciesData(2);

	DATA_STRUCTS::SpeciesData speciesA;
	DATA_STRUCTS::SpeciesData speciesB;

	speciesA.numParticles = numParticles;
	speciesA.spatialPerturbationMode = spatialPerturbationMode;
	speciesA.spatialPerturbationAmplitude = spatialPerturbationAmplitude;
	speciesA.driftVelocity = driftVelocity;
	speciesA.thermalVelocity = thermalVelocity;
	speciesA.plasmaFrequency = plasmaFrequency;
	speciesA.chargeMassRatio = chargeMassRatio;

	speciesB.numParticles = numParticles;
	speciesB.spatialPerturbationMode = spatialPerturbationMode;
	speciesB.spatialPerturbationAmplitude = spatialPerturbationAmplitude;
	speciesB.driftVelocity = - driftVelocity;
	speciesB.thermalVelocity = thermalVelocity;
	speciesB.plasmaFrequency = plasmaFrequency;
	speciesB.chargeMassRatio = chargeMassRatio;

	std::vector<double> APositions(speciesA.numParticles, 0);
	std::vector<double> AVelocities(speciesA.numParticles, 0);

	std::vector<double> BPositions(speciesB.numParticles, 0);
	std::vector<double> BVelocities(speciesB.numParticles, 0);

	speciesA.particlePositions = APositions;
	speciesA.particleXVelocities = AVelocities;

	speciesB.particlePositions = BPositions;
	speciesB.particleXVelocities = BVelocities;

	allSpeciesData[0] = speciesA;
	allSpeciesData[1] = speciesB;

	PIC_PLUS_PLUS::PICPlusPlus init(
		spatialLength,
		numTimeSteps,
		timeStepSize,
		numGrid,
		numSpecies,

		allSpeciesData);

	auto jsonResult = init.initialize();
	EXPECT_TRUE(jsonResult.has_value());

	EXPECT_EQ(jsonResult.value(), BuildExpectedJson());
}
