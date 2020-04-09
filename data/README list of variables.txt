List of variables

Column	Variable	Description

1	WM_ID		: unique grid cell ID
2	Centroid_X	: coordinates 
3	Centroid_Y	: coordinates
4	SR_M		: Species Richness mammals 
5	SR_A		: Species Richness amphibians 
6	SR_B		: Species Richness birds 
7	SR_Sum		: Species Richness total (sum mammals + amphibians + birds)

8	REGION_2016	: Regions as defined in 2016
9	REGION_2017	: Regions as defined in 2017

10	grid_bio1	: Annual mean temp (Chelsa)
11	grid_bio2	: Mean diurnal range (Chelsa)
12	grid_bio3	: Isothemality (Chelsa)
13	grid_bio4	: Temperature seasonality (Chelsa)
14	grid_bio5	: Max temp of warmest month (Chelsa)
15	grid_bio6	: Min temp of coldest month (Chelsa)
16	grid_bio7	: temp annual range (Chelsa)
17	grid_bio8	: Mean temp of wettest quarter (Chelsa)
18	grid_bio9	: Mean temp of driest quarter (Chelsa)
19	grid_bio10	: Mean temp of warmest quarter (Chelsa)
20	grid_bio11	: Mean temp of coldest quarter (Chelsa)
21	grid_bio12	: Annual precipitation (precip) (Chelsa)
22	grid_bio13	: Precip of wettest month (Chelsa)
23	grid_bio14	: Precip of driest month (Chelsa)
24	grid_bio15	: Precip seasonality (Chelsa)
25	grid_bio16	: Precip wettest quarter (Chelsa)
26	grid_bio17	: Precip driest quarter (Chelsa)
27	grid_bio18	: Precip of warmest quarter (Chelsa)
28	grid_bio19	: Precip of coldest quarter (Chelsa)
29	pet_ann		: A & T potential evapotranspiration


## GEOLOGICAL VARIABLES
30	Countpnts_t	: Number of geological datapoints on cooling age within a gridcell
31	Avg_t		: Cooling age (AFT) average (Ma)
32	Range_t		: Range of cooling age
33	Countpnts_SP	: Number of geological datapoints on erosion within a gridcell
34	Avg_TSP		: Total stream power index weight with precipitation data (Watts/m2)
35	Avg_USP		: Unit stream power index weight with precipitation data (Watts/m2)
36	Avg_SSP		: Shear stream power index weight with precipitation data (Watts/m2)
37	Avg_erlong	: Steady state long term erosion rate derived from Brandon's method (km/Ma)

38	Avg_HERM02	: Variable long-term erosion rate derived from Herman's method for interval between 2 and 0 Ma (km/Ma)
39	Avg_HERM24	: Variable long-term erosion rate derived from Herman's method for interval between 4 and 2 Ma (km/Ma)
40	Avg_HERM46	: Variable long-term erosion rate derived from Herman's method for interval between 6 and 4 Ma (km/Ma)
41	Avg_HERM68	: Variable long-term erosion rate derived from Herman's method for interval between 8 and 6 Ma (km/Ma)
42	Avg_HERM81	: Variable long-term erosion rate derived from Herman's method for interval between 10 and 8 Ma (km/Ma)
43	Avg_HERM10	: Variable long-term erosion rate derived from Herman's method for interval between 12 and 10 Ma (km/Ma)
44	Avg_HERMER	: Average variable long-term erosion rate derived from Herman's method for interval between 12 and 0 Ma (km/Ma)

## TOPOGRAPHIC VARIABLES
45	soil_var	: Soil variability within gridcell
46	dem_range	: Elevation range within gridcell
47 	Avg_Z		: Average elevation based on geological datapoints and calculated by Mauricio
48	Min_Z		: Minimal elevation based on geological datapoints and calculated by Mauricio
49	Max_Z		: Maximum elevation based on geological datapoints and calculated by Mauricio
50	Range_Z		: Range elevation based on geological datapoints and calculated by Mauricio
51	Rugg_avg	: Average ruggedness value based on Korner et al. (2011) Alp Bot.
52	Rugg_range	: Ruggedness range based on Korner et al (2011) Alp Bot.

53	Countpnts_cosmo	: Number of data values on cosmogenic isotope-derived millennial-scale erosion rates based on Willenbring.
54	E_calc		: Cosmogenic isotope-derived millennial-scale erosion rates (Willenbring)
55	dE_Calc		: SD (?) cosmogenic isotope-derived millennial-scale erosion rates (Willenbring)

56 	RELIEF_MIN     	: min value of range of elevation values (max-min) within a 2.5 km radius of a 90m resolution pixel
57 	RELIEF_MAX     	: max value of range of elevation values (max-min) within a 2.5 km radius of a 90m resolution pixel
58 	RELIEF_RNG     	: range of value of range of elevation values (max-min) within a 2.5 km radius of a 90m resolution pixel
59 	RELIEF_MN      	: mean range of elevation values (max-min) within a 2.5 km radius of a 90m resolution pixel
60 	RELIEF_STD     	: SD range of elevation values (max-min) within a 2.5 km radius of a 90m resolution pixel
61 	REALM   	: Realm regions

## LATER ADDED IN R-SCRIPT
62	HERMratio	: Avg_HERM02/Avg_HERM46
63	SR_All		: Species richness (same as SR_Sum)
