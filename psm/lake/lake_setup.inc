! Version 1.0
! Declarations ********************************************************
      real dz,dt,readdt,basedep,zo,z_screen,delta,rhowat,rhoice,rhosnow
      real Le,Lei,fusion,surf,fracmin,fraclim,qwtau,cpw_ice,condw,beta
      real snocrit,kv,dm,pi,tsurfp,mixing,trace,depth_begin,salty_begin
      real trace_i,trace_a,d_frac,tempi_a,hice_a,hsnow_a,xlat,xlon,gmt
      real b_area,runout_sum,alb_slush,alb_snow,fsol1,fsol2,sical0
      real sical1
      real salty_a,eta,fraci_a,temp_a,d_fraca,snowcut,o18air,f,ti_a
      real deutair,dexch,alphak,lamisw,lamilw,ca,cb,c70,c71,c72,c73
      real lamssw,lamslw,afrac1,afrac2,condi,conds,cpair,emis
      real o18snow_a,deutsnow_a,psum_a,surfarea_a,hsnow_ave,temp_ave
      real tsurf_ave,fice_ave,evap_ave,hice_ave,o18_ave,deut_ave
      real oblq,cdrn,raddeg,dpd,grav,sigma,rair,rvap,cvap
      real co18prec,co18run,cdeutprec,cdeutrun,xt,mix_ave

      integer max_dep,ix1,iy1,n_trace,i_area,lcount,nspin
      integer iin,iout,numpts,num_tra,ktau,iupto,ilake,jlake
      integer depth,depth_a,ktauwan,mixmax_a

      logical bndry_flag,iceflag,s_flag,o18flag,deutflag,ar_flag
      logical snow_flag_a,wb_flag,melt_flag_a

      real tempinit, deutinit, o18init !Ashling
!**********************************************************************
! Lake specific parameters ********************************************

      parameter (oblq = 23.4)           ! obliquity
      parameter (xlat = -6.30)          ! latitude (negative for South)
      parameter (xlon = 29.5)           ! longitude (negative for West)
      parameter (gmt = +3)              ! local time relative to gmt in hours
      parameter (max_dep = 999)         ! depth of lake at sill in meters 
      parameter (basedep = 733.)	    ! elevation of basin bottom in meters
      parameter (b_area = 23100000.)    ! area of catchment+lake in hectares 
      parameter (cdrn = 2.0e-3)         ! neutral drag coefficient 1.8 HAD 1.7GISS 1.2CCSM
      parameter (eta = 0.065)             ! shortwave extinction coefficient (1/m)
      parameter (f= 0.3)                 ! fraction of advected air over lake
      parameter (alb_slush = 0.4)       ! albedo of melting snow
      parameter (alb_snow = 0.7)        ! albedo of non-melting snow

      parameter (depth_begin = 570)       ! prescribed depth in meters
      parameter (salty_begin = 0.0)     ! prescribed salinity in ppt

      parameter (o18air = -14.0)         ! d18O of air above lake
      parameter (deutair = -96.)       ! dD of air above lake


      ! Lake area moved to bottom

      !real area(max_dep)
      !data area/570*3290000./           ! lake area in hectares by depth
     
      !integer, parameter :: hyplength=1500
      !real, dimension(hyplength):: area
      !open(unit=100, file='hypsograph.txt') 
      !read(100,*) area
      !close(unit=100) ! for orderly closing

      !Begin Ashling
      ! Removed from the model, for easier access
      character(38) :: datafile='ERA_INTERIM_climatology_Tanganyika.txt' ! the data file to open in FILE_OPEN subroutine
      parameter (tempinit= 23.0) ! temperature to initialize lake at in INIT_LAKE subroutine
      parameter (deutinit= 24.0)! dD to initialize lake at in INIT_LAKE subroutine
      parameter (o18init= 3.7)  ! d18O to initialize lake at in INIT_LAKE subroutine
      !End Ashling

! Simulation specific parameters **************************************

      parameter (nspin = 30)              ! number of years for spinup
      parameter (bndry_flag = .false.)    ! true for explict boundry layer computations; presently only for sigma coord climate models
      parameter (sigma = 0.9925561)       ! sigma level for boundary flag
      parameter (wb_flag = .true.)       ! true for variable lake depth
      parameter (iceflag = .false.)        ! true for variable ice cover
      parameter (s_flag = .false.)        ! true for variable salinity
      parameter (o18flag = .true.)       ! true for variable d18O 
      parameter (deutflag = .true.)      ! true for variable dD
      parameter (z_screen = 5.0)          ! height of met inputs

!**********************************************************************
! Other parameters DO NOT CHANGE without good reason for doing so******

      parameter (dz=1.0)                ! vertical layer thickness in m
      parameter (dt = 1.*60.*30.)       ! model time step in seconds
      parameter (zo=0.0004)             ! water surface roughness length
      parameter (n_trace=3)             ! model does 3 tracers (o18, D, sal)
      parameter (delta=5.67e-8)         ! s-b constant
      parameter (rhowat = 1000.)        ! density of water
      parameter (rhosnow = 330.)        ! density of snow
      parameter (rhoice = 917.)         ! density of ice
      parameter (snowcut = 2.2)         ! temp at which precip is snow
      parameter (cpair = 1.00464e3)     ! specific heat capacity dry air
      parameter (rair = 287.05)         ! specific gas constant dry air
      parameter (rvap = 461.495)        ! specific gas constant water vapor
      parameter (cvap = 1.81e3)         ! specific heat capacity water vapor
      parameter (Le = 2.45e6)           ! latent heat of vaporization water
      parameter (Lei = 2.5e6)           ! latent heat of vaporization ice
      parameter (ca = 6.11e2)           ! constant for Teten's formula
      parameter (cb = 273.16)           ! constant for Teten's formula
      parameter (c72 = 17.269)          ! constant for Teten over water
      parameter (c73 = 35.86)           ! constant for Teten over water
      parameter (c70 = 21.874)          ! constant for Teten over ice
      parameter (c71 = 7.66)            ! constant for Teten over ice
      parameter (emis = 0.97)           ! longwave emmisivity
      parameter (fusion=3.34e5)         ! latent heat of fusion
      parameter (surf = 1.0)            ! surface thickness over which the flux exchange occurs in water
      parameter (fracmin= 0.01)         ! min ice thick in meters
      parameter (fraclim = 0.02)        ! min ice fraction
      parameter (qwtau = 86400.)        ! D. Pollard sub-ice time constant
      parameter (cpw_ice = 4200.)       ! heat capacity of ice
      parameter (condw = .58)           ! conduction coefficinet for water
      parameter (beta=0.4)              ! frac of solar rad abs in water surface layer
      parameter (fsol1=0.5)             ! sw weighting for surface albedo
      parameter (fsol2=0.5)             ! lw weighting for surface albedo
      parameter (sical0=0.6)            ! sw albedo for sea ice
      parameter (sical1=0.4)            ! lw albedo for sea ice
      parameter (kv=0.4)                ! vonkarman constant
      parameter (dm=1.38889E-07)        ! thermal molecular diffusivity of water
      parameter (pi=3.141592654)        !
      parameter (snocrit = 0.05)        ! for albedo, in m
      parameter (dexch  =1.0)           ! fixed depth (m) over which isotope exchange occurs
      parameter (alphak = 0.994)        ! fractionation coeff
      parameter (lamisw = 3.0)          ! extinction coef for sw thru ice
      parameter (lamilw = 20.)          ! extinction coef for lw thru ice
      parameter (lamssw = 3.0)          ! extinction coef for sw thru snow
      parameter (lamslw = 20.)          ! extinction coef for lw thru snow
      parameter (afrac1 = 0.7)          ! fraction of light in visible band for ice model
      parameter (afrac2 = 0.3)          ! fraction of light in infrared band for ice model
      parameter (condi = 2.3)           ! thermal conductivity of ice
      parameter (conds = 0.31)          ! thermal conductivity of snow
      parameter (raddeg = 0.0174533)    ! radians per degree
      parameter (dpd = 360./365.)       ! degrees per day
      parameter (grav = 9.80616)        ! gravity

!**********************************************************************

! Common definitions
      common /ldata/ lcount, iin, iout, numpts , num_tra, ktau, iupto
      common mixing, o18snow_a,deutsnow_a
      common psum_a, depth_a
      common trace(max_dep,n_trace), trace_i(max_dep,n_trace)
      common trace_a (max_dep, n_trace),surfarea_a
      common d_fraca, tempi_a, hice_a, hsnow_a
      common salty_a, mixmax_a, fraci_a 
      common temp_a (max_dep), snow_flag_a, melt_flag_a
      common ktauwan (12),ti_a(max_dep)
      common mix_ave, tsurf_ave, fice_ave
      common evap_ave, hice_ave, hsnow_ave, o18_ave
      common deut_ave, runout_sum, temp_ave(max_dep)

!**********************************************************************
      !Ashling
      real, dimension(999) :: area
      data area/3290000,3246809,3235562,  &
      3227075,3219947,3213619,  &
      3207911,3202743,3197948,  &
      3193338,3188901,3184763,  &
      3180759,3176884,3173213,  &
      3169402,3165770,3162305,  &
      3158823,3155504,3151633,  &
      3147965,3144780,3141458,  &
      3138169,3134920,3131680,  &
      3128600,3125470,3122458,  &
      3119361,3116264,3113250,  &
      3110231,3107155,3104048,  &
      3100808,3097635,3094222,  &
      3090791,3086891,3083211,  &
      3079875,3076551,3073218,  &
      3069895,3066548,3063255,  &
      3059773,3056249,3051980,  &
      3047720,3044308,3041014,  &
      3037787,3034670,3031583,  &
      3028538,3025529,3022673,  &
      3019459,3016373,3013529,  &
      3010694,3007868,3005091,  &
      3002249,2999488,2996796,  &
      2993993,2991234,2988320,  &
      2985593,2982797,2980008,  &
      2977115,2974295,2971594,  &
      2968793,2966023,2962862,  &
      2959644,2956914,2954219,  &
      2951607,2948952,2946271,  &
      2943708,2941094,2938531,  &
      2935822,2933124,2930490,  &
      2927827,2925134,2922339,  &
      2919593,2916844,2914098,  &
      2911225,2907195,2903030,  &
      2900256,2897525,2894890,  &
      2892295,2889543,2886990,  &
      2884535,2882136,2879505,  &
      2876973,2874619,2872208,  &
      2869853,2867386,2864927,  &
      2862476,2860053,2857578,  &
      2854812,2852016,2849537,  &
      2847005,2844493,2841816,  &
      2839093,2836582,2834102,  &
      2831613,2828976,2826403,  &
      2823847,2821294,2818706,  &
      2815978,2813257,2810534,  &
      2807774,2804988,2802002,  &
      2798948,2795989,2792856,  &
      2789736,2786290,2782692,  &
      2779057,2775216,2771149,  &
      2765965,2761472,2758121,  &
      2754814,2751678,2748429,  &
      2745299,2742340,2739416,  &
      2736525,2733312,2730198,  &
      2727338,2724493,2721624,  &
      2718664,2715843,2713137,  &
      2710451,2707827,2705074,  &
      2702320,2699609,2697023,  &
      2694307,2691420,2688536,  &
      2685678,2682873,2679982,  &
      2676746,2673719,2670837,  &
      2667905,2664960,2662007,  &
      2659033,2656164,2653176,  &
      2650270,2647133,2643974,  &
      2640922,2637719,2634337,  &
      2630823,2627209,2623465,  &
      2619478,2614907,2608009,  &
      2602009,2597819,2593879,  &
      2590066,2586545,2583280,  &
      2580273,2577382,2574587,  &
      2571832,2569042,2566377,  &
      2563864,2561260,2558705,  &
      2556098,2553637,2551249,  &
      2548788,2546291,2543735,  &
      2541396,2539040,2536659,  &
      2534186,2531760,2529339,  &
      2526947,2524612,2522182,  &
      2519756,2517388,2514999,  &
      2512674,2510209,2507727,  &
      2505308,2502828,2500382,  &
      2497649,2494883,2492313,  &
      2489665,2486958,2484041,  &
      2481093,2478105,2474962,  &
      2471819,2467415,2462943,  &
      2459732,2456695,2453592,  &
      2450498,2447700,2445039,  &
      2442352,2439824,2437172,  &
      2434522,2432023,2429540,  &
      2427132,2424632,2422221,  &
      2419861,2417472,2415107,  &
      2412745,2410443,2408186,  &
      2405846,2403590,2401177,  &
      2398846,2396506,2394200,  &
      2391857,2389390,2387018,  &
      2384659,2382360,2379964,  &
      2377405,2374908,2372521,  &
      2370143,2367645,2365028,  &
      2362451,2359744,2357081,  &
      2354298,2351306,2348188,  &
      2344983,2341586,2337713,  &
      2330287,2321769,2317117,  &
      2313063,2309294,2305612,  &
      2302077,2298786,2295646,  &
      2292527,2289409,2286354,  &
      2283296,2280271,2277237,  &
      2274005,2270728,2267382,  &
      2264142,2261082,2258157,  &
      2255125,2252243,2249313,  &
      2246576,2243680,2240812,  &
      2238121,2235423,2232700,  &
      2230025,2227365,2224766,  &
      2222130,2219581,2217080,  &
      2214638,2212179,2209764,  &
      2207197,2204543,2201823,  &
      2199108,2196277,2193484,  &
      2190538,2187464,2184000,  &
      2180302,2176056,2169272,  &
      2162757,2158476,2154595,  &
      2150874,2147176,2143568,  &
      2140203,2136938,2133868,  &
      2130647,2127582,2124689,  &
      2121926,2119367,2116644,  &
      2113962,2111353,2108881,  &
      2106290,2103725,2101225,  &
      2098763,2096362,2093953,  &
      2091552,2089064,2086579,  &
      2084050,2081451,2078706,  &
      2075998,2073360,2070728,  &
      2068182,2065508,2062881,  &
      2060209,2057522,2054805,  &
      2052041,2049191,2046261,  &
      2043315,2040284,2037282,  &
      2034073,2030728,2027326,  &
      2023610,2018527,2013427,  &
      2009730,2006196,2002742,  &
      1999303,1996171,1993094,  &
      1990087,1987022,1983867,  &
      1980823,1977813,1974825,  &
      1971978,1969156,1966311,  &
      1963510,1960692,1957944,  &
      1955109,1952284,1949778,  &
      1947248,1944839,1942361,  &
      1939985,1937652,1935338,  &
      1932974,1930638,1928272,  &
      1925972,1923600,1921180,  &
      1918811,1916356,1913810,  &
      1911289,1908772,1906012,  &
      1903179,1900462,1897769,  &
      1895019,1892302,1889430,  &
      1886375,1883332,1880088,  &
      1875679,1871143,1867741,  &
      1864556,1861342,1858223,  &
      1855249,1852260,1849273,  &
      1846133,1842547,1838867,  &
      1835994,1833270,1830820,  &
      1828348,1825904,1823459,  &
      1821054,1818780,1816459,  &
      1814141,1811771,1809438,  &
      1807055,1804692,1802311,  &
      1799838,1797262,1794510,  &
      1791498,1788159,1785255,  &
      1782536,1779865,1777296,  &
      1774691,1772161,1769620,  &
      1767073,1764480,1761892,  &
      1759190,1756542,1753789,  &
      1750978,1747992,1744926,  &
      1741685,1738000,1732322,  &
      1727192,1723322,1719961,  &
      1716686,1713521,1710455,  &
      1707491,1704534,1701613,  &
      1698803,1695987,1693277,  &
      1690485,1687755,1684959,  &
      1682216,1679383,1676513,  &
      1673555,1670360,1667182,  &
      1664126,1661105,1658078,  &
      1655084,1652092,1649005,  &
      1645992,1642821,1639675,  &
      1636485,1633204,1629816,  &
      1626396,1623036,1619621,  &
      1616134,1612533,1608917,  &
      1604880,1600935,1596931,  &
      1592829,1588553,1584233,  &
      1579895,1575396,1570579,  &
      1565211,1556220,1547382,  &
      1540840,1534495,1528629,  &
      1523216,1518725,1514838,  &
      1511039,1507454,1503784,  &
      1500291,1497112,1494025,  &
      1491136,1488261,1485453,  &
      1482732,1480246,1477875,  &
      1475436,1473059,1470752,  &
      1468367,1466057,1463819,  &
      1461495,1459214,1456920,  &
      1454553,1452052,1449372,  &
      1447052,1444756,1442573,  &
      1440399,1438328,1436213,  &
      1434106,1432025,1429959,  &
      1427870,1425828,1423680,  &
      1421504,1419213,1416990,  &
      1414525,1411973,1409286,  &
      1405701,1402026,1399330,  &
      1396674,1394128,1391710,  &
      1389321,1387016,1384909,  &
      1382736,1380756,1378732,  &
      1376685,1374670,1372611,  &
      1370319,1367976,1365653,  &
      1363400,1360956,1358690,  &
      1356537,1354582,1352729,  &
      1350908,1349196,1347435,  &
      1345722,1344017,1342356,  &
      1340690,1339010,1337407,  &
      1335768,1334122,1332419,  &
      1330795,1329167,1327439,  &
      1325709,1323938,1322109,  &
      1320226,1318405,1316708,  &
      1314957,1313121,1311271,  &
      1309473,1307559,1305320,  &
      1302913,1301081,1299394,  &
      1297723,1295993,1294443,  &
      1292944,1291351,1289832,  &
      1288350,1286873,1285401,  &
      1283926,1282460,1281026,  &
      1279556,1278106,1276721,  &
      1275290,1273846,1272453,  &
      1271046,1269600,1268158,  &
      1266767,1265394,1263954,  &
      1262571,1261186,1259761,  &
      1258293,1256882,1255502,  &
      1254032,1252587,1251118,  &
      1249723,1248160,1246703,  &
      1245245,1243634,1242149,  &
      1240558,1238938,1237252,  &
      1235666,1233875,1232051,  &
      1230292,1227605,1225165,  &
      1223475,1221674,1219946,  &
      1218282,1216574,1214924,  &
      1213363,1211744,1210206,  &
      1208671,1207135,1205619,  &
      1204101,1202558,1201113,  &
      1199601,1198057,1196565,  &
      1195066,1193563,1192002,  &
      1190502,1188989,1187440,  &
      1185916,1184359,1182805,  &
      1181254,1179607,1178026,  &
      1176438,1174812,1173201,  &
      1171446,1169852,1168157,  &
      1166459,1164726,1163008,  &
      1161240,1159406,1157471,  &
      1155558,1153504,1151404,  &
      1149334,1146981,1144677,  &
      1141569,1138529,1136338,  &
      1133981,1131991,1129880,  &
      1127858,1125903,1123924,  &
      1122039,1120178,1118346,  &
      1116492,1114685,1112833,  &
      1111024,1109151,1107362,  &
      1105552,1103738,1101913,  &
      1100074,1098272,1096426,  &
      1094506,1092667,1090739,  &
      1088939,1087028,1085158,  &
      1083308,1081432,1079514,  &
      1077559,1075657,1073721,  &
      1071712,1069745,1067788,  &
      1065680,1063721,1061625,  &
      1059492,1057416,1055217,  &
      1053018,1050756,1048438,  &
      1046019,1043587,1040492,  &
      1037133,1034536,1031982,  &
      1029450,1026945,1024522,  &
      1022102,1019697,1017335,  &
      1014985,1012641,1010253,  &
      1007892,1005412,1002815,  &
      1000423,998109,995884,  &
      993705,991509,989445,  &
      987359,985307,983307,  &
      981316,979352,977579,  &
      976091,974587,973148,  &
      971682,970201,968740,  &
      967319,965858,964431,  &
      963090,961669,960262,  &
      958879,957493,956110,  &
      954653,953277,951890,  &
      950704,949531,948368,  &
      947163,945836,944400,  &
      943123,941817,940541,  &
      939325,938090,936929,  &
      935864,934755,933698,  &
      932679,931644,930610,  &
      929605,928601,927590,  &
      926607,925570,924599,  &
      923593,922584,921614,  &
      920593,919545,918549,  &
      917552,916502,915502,  &
      914476,913444,912403,  &
      911326,910224,909131,  &
      908073,907004,905888,  &
      904697,903544,902394,  &
      901124,899696,898292,  &
      896967,895786,894433,  &
      892921,891473,889954,  &
      887203,884429,882897,  &
      881498,880125,878838,  &
      877694,876537,875325,  &
      874286,873297,872265,  &
      871315,870364,869390,  &
      868426,867474,866485,  &
      865575,864630,863689,  &
      862718,861815,860843,  &
      859932,858968,857957,  &
      856985,856040,855072,  &
      854128,853098,852103,  &
      851073,850028,848962,  &
      847915,846838,845740,  &
      844585,843485,842302,  &
      841117,839878,838655,  &
      837315,835908,834481,  &
      832885,831156,828871,  &
      826733,825110,823530,  &
      822064,820675,819316,  &
      818038,816738,815477,  &
      814251,812961,811728,  &
      810515,809294,808126,  &
      806987,805768,804608,  &
      803416,802302,801183,  &
      800095,798953,797830,  &
      796764,795613,794478,  &
      793375,792284,791213,  &
      790079,788999,787901,  &
      786854,785689,784578,  &
      783426,782259,781059,  &
      779902,778748,777497,  &
      776261,775027,773730,  &
      772186,770693,769116/!,  &
      !767464/
      