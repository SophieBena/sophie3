//
//  ZhangAndDuanCorrectionTerms.h
//  ThermoFit
//
//  Created by Mark Ghiorso on 4/15/16.
//  Copyright © 2016 Mark Ghiorso. All rights reserved.
//

#ifndef ZhangAndDuanCorrectionTerms_h
#define ZhangAndDuanCorrectionTerms_h

typedef struct _zAndDcorrections {
    double t;
    double g;
    double h;
    double s;
    double cp;
    double dcpdt;
} ZandDcorrections;

const static ZandDcorrections zAndDcorrections[] = {
    { 270, 45039.85053, 163354.1531, 438.2011291, -1486.30666, 17.92271899 },
    { 280, 40919.27332, 149343.9978, 387.2311428, -1319.863346, 15.43848334 },
    { 290, 37271.18718, 136880.9981, 343.482104, -1176.211447, 13.35210746 },
    { 300, 34029.611, 125756.01, 305.75467, -1051.7193, 11.592283 },
    { 310, 31139.313, 115792.62, 273.0752, -943.43561, 10.10518 },
    { 320, 28553.945, 106841.62, 244.64898, -848.87035, 8.8417461 },
    { 330, 26234.341, 98776.334, 219.82422, -765.98133, 7.7642901 },
    { 340, 24147.258, 91488.771, 198.06327, -693.06794, 6.8419677 },
    { 350, 22264.357, 84886.483, 178.92036, -628.70984, 6.049418 },
    { 360, 20561.369, 78890.045, 162.0241, -571.71678, 5.3658192 },
    { 370, 19017.428, 73430.952, 147.06358, -521.08773, 4.7740205 },
    { 380, 17614.521, 68449.908, 133.77733, -475.97773, 4.2598701 },
    { 390, 16337.04, 63895.398, 121.94451, -435.67093, 3.8116101 },
    { 400, 15171.413, 59722.517, 111.37776, -399.55869, 3.4194793 },
    { 410, 14105.799, 55891.982, 101.91752, -367.12165, 3.0753259 },
    { 420, 13129.832, 52369.323, 93.42736, -337.91513, 2.772325 },
    { 430, 12234.41, 49124.191, 85.790187, -311.5571, 2.5047538 },
    { 440, 11411.522, 46129.788, 78.90515, -287.71828, 2.2677528 },
    { 450, 10654.092, 43362.381, 72.685087, -266.1139, 2.0572376 },
    { 460, 9955.8569, 40800.889, 67.054417, -246.49696, 1.8697293 },
    { 470, 9311.2591, 38426.536, 61.947398, -228.65252, 1.7022843 },
    { 480, 8715.3554, 36222.558, 57.306672, -212.39295, 1.5523746 },
    { 490, 8163.7394, 34173.944, 53.082051, -197.55398, 1.4178481 },
    { 500, 7652.4748, 32267.226, 49.229502, -183.99128, 1.2968587 },
    { 510, 7178.0387, 30490.29, 45.710296, -171.57766, 1.1876457 },
    { 520, 6737.272, 28832.219, 42.490283, -160.20059, 1.0891408 },
    { 530, 6327.3367, 27283.157, 39.539285, -149.76013, 1.0000356 },
    { 540, 5945.679, 25834.193, 36.830582, -140.16711, 0.9193478 },
    { 550, 5589.9967, 24477.258, 34.340475, -131.34153, 0.846215 },
    { 560, 5258.2112, 23205.046, 32.047919, -123.2111, 0.7799196 },
    { 570, 4948.4424, 22010.939, 29.934204, -115.70995, 0.7198632 },
    { 580, 4658.9866, 20888.952, 27.982699, -108.77733, 0.6655649 },
    { 590, 4388.2972, 19833.687, 26.178627, -102.35631, 0.6166528 },
    { 600, 4134.9665, 18840.302, 24.508894, -96.392376, 0.5728691 },
    { 610, 3897.7102, 17904.495, 22.961942, -90.831803, 0.5340626 },
    { 620, 3675.3522, 17022.505, 21.527665, -85.619549, 0.5001807 },
    { 630, 3466.7324, 16190.873, 20.197049, -80.759613, 0.4555729 },
    { 640, 3271.0205, 15405.433, 18.960019, -76.398056, 0.413445 },
    { 650, 3087.2512, 14662.05, 17.807383, -72.308197, 0.391237 },
    { 660, 2914.6103, 13958.72, 16.733499, -68.386476, 0.369032 },
    { 670, 2752.3353, 13293.689, 15.733364, -64.654187, 0.3447 },
    { 680, 2599.7147, 12664.9, 14.801743, -61.141838, 0.320853 },
    { 690, 2425.3585, 12126.76, 14.060002, -55.205467, 0.299928 },
    { 700, 2288.6708, 11589.025, 13.28622, -52.363481, 0.282965 },
    { 710, 2159.4661, 11079.091, 12.562853, -49.641338, 0.269871 },
    { 720, 2037.2555, 10595.872, 11.886968, -47.017116, 0.259967 },
    { 730, 1921.5769, 10138.486, 11.25604, -44.471907, 0.252337 },
    { 740, 1811.9922, 9706.2153, 10.667869, -41.992614, 0.245864 },
    { 750, 1708.0837, 9298.4262, 10.120457, -39.576127, 0.239112 },
    { 760, 1609.4537, 8914.4454, 9.611831, -37.234529, 0.230222 },
    { 770, 1515.725, 8553.3814, 9.139814, -34.999948, 0.2170172 },
    { 780, 1426.5442, 8213.9122, 8.701754, -32.926219, 0.197453 },
    { 790, 1341.5878, 7894.0871, 8.294302, -31.083562, 0.1704893 },
    { 800, 1260.5698, 7591.2297, 7.913325, -29.543588, 0.1370493 },
    { 810, 1183.2489, 7302.0335, 7.554055, -28.356667, 0.1004524 },
    { 820, 1109.4329, 7022.8885, 7.211531, -27.53029, 0.06572863 },
    { 830, 1038.9778, 6750.3691, 6.881195, -27.020039, 0.0378869 },
    { 840, 971.78056, 6481.7207, 6.559452, -26.739252, 0.0201114 },
    { 850, 907.76791, 6215.1715, 6.244004, -26.582479, 0.012818 },
    { 860, 846.88296, 5949.9905, 5.933846, -26.450714, 0.014598 },
    { 870, 789.07326, 5686.3292, 5.62903, -26.268364, 0.022487 },
    { 880, 734.2821, 5424.9497, 5.330304, -25.988983, 0.033602 },
    { 890, 682.44324, 5166.942, 5.038762, -25.592517, 0.045609 },
    { 900, 633.479, 4913.4924, 4.755571, -25.078549, 0.056923 },
    { 910, 587.30046, 4665.7236, 4.481784, -24.459007, 0.066639 },
    { 920, 543.80911, 4424.6035, 4.218255, -23.752141, 0.074368 },
    { 930, 502.8991, 4190.9041, 3.965597, -22.978269, 0.0800603 },
    { 940, 464.45965, 3965.1955, 3.724187, -22.157128, 0.0838613 },
    { 950, 428.37729, 3747.8597, 3.494192, -21.306448, 0.0860158 },
    { 960, 394.53779, 3539.1143, 3.275601, -20.441341, 0.0867958 },
    { 970, 362.82777, 3339.0394, 3.068259, -19.574185, 0.0864701 },
    { 980, 333.13596, 3147.6044, 2.871907, -18.714796, 0.0852809 },
    { 990, 305.35414, 2964.6922, 2.6862, -17.870725, 0.0834389 },
    { 1000, 279.37779, 2790.1199, 2.510742, -17.047597, 0.0811196 },
    { 1010, 255.10658, 2623.6569, 2.345099, -16.249444, 0.0784656 },
    { 1020, 232.4446, 2465.0385, 2.188818, -15.47902, 0.0755912 },
    { 1030, 211.30057, 2313.9782, 2.041435, -14.738058, 0.0725864 },
    { 1040, 191.58779, 2170.179, 1.902488, -14.027497, 0.0695209 },
    { 1050, 173.22419, 2033.3258, 1.771525, -13.347662, 0.0664483 },
    { 1060, 156.13215, 1903.1207, 1.648102, -12.698419, 0.0634084 },
    { 1070, 140.23843, 1779.257, 1.531793, -12.079288, 0.0604301 },
    { 1080, 125.47394, 1661.437, 1.422188, -11.48954, 0.0575346 },
    { 1090, 111.77362, 1549.3712, 1.318896, -10.928271, 0.0547363 },
    { 1100, 99.076198, 1442.78, 1.221549, -10.394457, 0.0520446 },
    { 1110, 87.324008, 1341.3942, 1.129793, -9.887003, 0.0494652 },
    { 1120, 76.462805, 1244.9559, 1.043298, -9.404771, 0.0470005 },
    { 1130, 66.441569, 1153.2186, 0.961749, -8.94661, 0.0446509 },
    { 1140, 57.212304, 1065.9473, 0.884856, -8.511372, 0.0424154 },
    { 1150, 48.729875, 982.91848, 0.812337, -8.09793, 0.04029147 },
    { 1160, 40.951824, 903.91971, 0.743938, -7.705183, 0.03827585 },
    { 1170, 33.838209, 828.7494, 0.679411, -7.332065, 0.03636472 },
    { 1180, 27.351451, 757.2164, 0.61853, -6.977553, 0.03455395 },
    { 1190, 21.456183, 689.1396, 0.561079, -6.640665, 0.03283914 },
    { 1200, 16.119117, 624.34747, 0.506857, -6.320466, 0.03121563 },
    { 1210, 11.308907, 562.67763, 0.455676, -6.016062, 0.02967904 },
    { 1220, 6.996036, 503.97639, 0.407361, -5.726609, 0.02822489 },
    { 1230, 3.152692, 448.09829, 0.361744, -5.451304, 0.02684882 },
    { 1240, -0.24733, 394.90568, 0.318672, -5.189387, 0.02554652 },
    { 1250, -3.228737, 344.26831, 0.277997, -4.940141, 0.02431395 },
    { 1260, -5.814822, 296.06289, 0.239586, -4.702888, 0.02314733 },
    { 1270, -8.027559, 250.17271, 0.203308, -4.476988, 0.02204285 },
    { 1280, -9.88768, 206.48731, 0.169043, -4.261836, 0.02099694 },
    { 1290, -11.414754, 164.9027, 0.13668, -4.056864, 0.02000638 },
    { 1300, -12.627258, 125.3179, 0.106112, -3.861534, 0.019067832 },
    { 1310, -13.54264, 87.640922, 0.07724, -3.675342, 0.018178488 },
    { 1320, -14.177385, 51.782183, 0.049969, -3.497811, 0.017335461 },
    { 1330, -14.547075, 17.657352, 0.024214, -3.328488, 0.016536028 },
    { 1340, -14.666435, -14.813535, -0.000109, -3.166953, 0.015777733 },
    { 1350, -14.549392, -45.706326, -0.023079, -3.012804, 0.015058208 },
    { 1360, -14.209121, -75.092985, -0.044768, -2.865666, 0.01437538 },
    { 1370, -13.658082, -103.04182, -0.065244, -2.725182, 0.01372699 },
    { 1380, -12.908073, -129.61769, -0.084572, -2.591018, 0.013111092 },
    { 1390, -11.970259, -154.88219, -0.102814, -2.462857, 0.01252598 },
    { 1400, -10.855215, -178.89385, -0.120027, -2.340401, 0.01196996 },
    { 1410, -9.572955, -201.70829, -0.136266, -2.223368, 0.01144119 },
    { 1420, -8.132969, -223.37839, -0.151581, -2.111491, 0.01093835 },
    { 1430, -6.544248, -243.95445, -0.166022, -2.004519, 0.01045995 },
    { 1440, -4.815313, -263.48433, -0.179632, -1.902215, 0.01000466 },
    { 1450, -2.954245, -282.01355, -0.192455, -1.804353, 0.00957121 },
    { 1460, -0.968705, -299.58549, -0.204532, -1.710721, 0.00915842 },
    { 1470, 1.134037, -316.24142, -0.215902, -1.621119, 0.0087652 },
    { 1480, 3.34709, -332.0267, -0.226599, -1.535356, 0.00839048 },
    { 1490, 5.663911, -346.9673, -0.236661, -1.453252, 0.00803328 },
    { 1500, 8.078291, -361.09733, -0.246117, -1.374635, 0.00769269 },
    { 1510, 10.584339, -374.46452, -0.254999, -1.299346, 0.00736782 },
    { 1520, 13.17646, -387.09481, -0.263337, -1.227229, 0.00705786 },
    { 1530, 15.849342, -399.0192, -0.271156, -1.158141, 0.00676204 },
    { 1540, 18.597939, -410.26727, -0.278484, -1.091944, 0.00647963 },
    { 1550, 21.417463, -420.86727, -0.285345, -1.028506, 0.00620994 },
    { 1560, 24.303361, -430.84617, -0.291762, -0.967704, 0.00595233 },
    { 1570, 27.251309, -440.22975, -0.297759, -0.909421, 0.0057062 },
    { 1580, 30.257198, -449.04261, -0.303355, -0.853544, 0.00547096 },
    { 1590, 33.317127, -457.30829, -0.30857, -0.799967, 0.00524608 },
    { 1600, 36.427382, -465.04929, -0.313423, -0.74859, 0.00503103 },
    { 1610, 39.584439, -472.2871, -0.317932, -0.699316, 0.00482534 },
    { 1620, 42.784946, -479.0423, -0.322116, -0.652053, 0.00462855 },
    { 1630, 46.025716, -485.33458, -0.325988, -0.606716, 0.00444023 },
    { 1640, 49.30372, -491.18277, -0.329565, -0.563222, 0.004259965 },
    { 1650, 52.616079, -496.60489, -0.332862, -0.521492, 0.004087377 },
    { 1660, 55.960054, -501.61822, -0.335891, -0.48145, 0.003922098 },
    { 1670, 59.333043, -506.23928, -0.338666, -0.443026, 0.003763779 },
    { 1680, 62.732571, -510.48391, -0.3412, -0.406152, 0.0036121 },
    { 1690, 66.156284, -514.36727, -0.343505, -0.370763, 0.003466744 },
    { 1700, 69.601945, -517.90391, -0.345592, -0.336797, 0.003327412 },
    { 1710, 73.067426, -521.10776, -0.34747, -0.304195, 0.003193832 },
    { 1720, 76.550705, -523.99218, -0.349153, -0.272902, 0.003065745 },
    { 1730, 80.049858, -526.56998, -0.350648, -0.242863, 0.002942895 },
    { 1740, 83.563057, -528.85346, -0.351964, -0.214027, 0.002825042 },
    { 1750, 87.088564, -530.85438, -0.353111, -0.186346, 0.00271196 },
    { 1760, 90.624725, -532.58408, -0.354096, -0.159773, 0.00260345 },
    { 1770, 94.169972, -534.05339, -0.354928, -0.134263, 0.00249929 },
    { 1780, 97.722809, -535.27273, -0.355615, -0.109773, 0.00239929 },
    { 1790, 101.28182, -536.25212, -0.356164, -0.086264, 0.00230327 },
    { 1800, 104.84565, -537.00114, -0.356582, -0.063695, 0.00221107 },
    { 1810, 108.41303, -537.52903, -0.356874, -0.04203, 0.00212249 },
    { 1820, 111.98274, -537.84464, -0.357048, -0.021234, 0.0020374 },
    { 1830, 115.55362, -537.95649, -0.357109, -0.001271, 0.00195563 },
    { 1840, 119.12457, -537.87274, -0.357064, 0.01789, 0.00187706 },
    { 1850, 122.69456, -537.60126, -0.356917, 0.03628, 0.00180154 },
    { 1860, 126.26258, -537.14961, -0.356673, 0.05393, 0.00172894 },
    { 1870, 129.82772, -536.52503, -0.356338, 0.070868, 0.00165914 },
    { 1880, 133.38906, -535.73452, -0.355917, 0.087121, 0.00159202 },
    { 1890, 136.94578, -534.78479, -0.355413, 0.102717, 0.00152746 },
    { 1900, 140.49706, -533.68229, -0.354831, 0.117679, 0.00146538 },
    { 1910, 144.04216, -532.43324, -0.354176, 0.132032, 0.00140564 },
    { 1920, 147.58034, -531.04359, -0.35345, 0.1458, 0.00134818 },
    { 1930, 151.11094, -529.51912, -0.352658, 0.159003, 0.00129288 },
    { 1940, 154.63329, -527.86533, -0.351803, 0.171664, 0.00123967 },
    { 1950, 158.14681, -526.08757, -0.35089, 0.183804, 0.00118846 },
    { 1960, 161.65089, -524.1994, -0.349919, 0.19544, 0.00113915 },
    { 1970, 165.14502, -522.1838, -0.348896, 0.206592, 0.00109169 },
    { 1980, 168.62865, -520.06, -0.347822, 0.217279, 0.00104599 },
    { 1990, 172.10132, -517.83628, -0.346702, 0.227518, 0.00100199 },
    { 2000, 175.56255, -515.51172, -0.345537, 0.237325, 0.00095961 }
};

const static NSInteger nZandDcorrections = (sizeof(zAndDcorrections)/sizeof(struct _zAndDcorrections));

#endif /* ZhangAndDuanCorrectionTerms_h */