from scipy.stats import norm 
import numpy as np 
from pines_analysis_toolkit.astrometry import faherty_2016_spt_teff_relation
import astropy.constants as const 

def pines_J_band_noise_model(x, nexp, exp_time):
    '''Given a magnitude, nexp, and exp_time, calculate the expected stddev of photometric measurements using the noise model from Tamburo et al. (2022).
    x: J-band magnitude
    nexp: number of exposures to bin 
    exp_time: exposure time of individual exposures in seconds
    '''

    seeing = 2.5 
    pscale = 0.579
    aperture_radius = 2.0 
    dark_current = 22.65 /exp_time
    read_noise = 19
    central_blockage = 0.05
    net_detection = 8.21006209e-02
    L = 7.87029106e-03
    k = 1.00000000e+00
    x0 =  1.44492793e+01

    D_tel = 1.83 #diameter of telescope in meters
    tel_A = (1 - central_blockage) * np.pi * (D_tel * 100. / 2.) ** 2.

    phi_sky = 17.9*8.21 #e- pix^-1 s^-1
    phi = 3.149e5  # [photons cm^-2 s^-1], these are in the 2mass photometric system
    sigma_test = aperture_radius * seeing / 2.355
    area_under_gaussian = norm(0, seeing / 2.355).cdf(sigma_test) - norm(0, seeing / 2.355).cdf(-sigma_test)
    radius = aperture_radius * seeing / pscale  # Radius of this aperture in pixels
    npix = np.pi * radius ** 2.  # Number of pixels contained within the aperture
    
    phi_star = phi * tel_A * 10 ** (-x / 2.5) * net_detection * area_under_gaussian  # [e- s^-1]
    signal = phi_star * exp_time  # [e-]
    noise = np.sqrt(signal + npix * (phi_sky * exp_time + dark_current * exp_time + read_noise ** 2.))
    snr = (signal / noise)
    #Add a logistic noise floor
    output = (1/snr)/np.sqrt(nexp) + L*(1-1/(1+np.exp(-k*(x-x0))))
    return output

def baraffe_models():
    ages = np.array([5.00e-04, 1.00e-03, 2.00e-03, 3.00e-03, 4.00e-03, 5.00e-03,
       8.00e-03, 1.00e-02, 1.50e-02, 2.00e-02, 2.50e-02, 3.00e-02,
       4.00e-02, 5.00e-02, 8.00e-02, 1.00e-01, 1.20e-01, 2.00e-01,
       3.00e-01, 4.00e-01, 5.00e-01, 6.25e-01, 8.00e-01, 1.00e+00,
       2.00e+00, 3.00e+00, 4.00e+00, 5.00e+00, 8.00e+00, 1.00e+01])
    
    t_effs = np.array([list([2387.0, 2514.0, 2594.0, 2668.0, 2733.0, 2762.0, 2808.0, 2835.0, 2832.0, 2836.0, 2837.0, 2864.0, 2925.0, 2945.0, 3006.0, 3076.0, 3139.0, 3220.0, 3460.0, 3670.0, 3848.0, 3987.0, 4110.0, 4221.0, 4315.0, 4397.0, 4471.0, 4537.0, 4596.0, 4647.0]),
                       list([2345.0, 2504.0, 2598.0, 2710.0, 2779.0, 2824.0, 2864.0, 2897.0, 2896.0, 2897.0, 2898.0, 2908.0, 2936.0, 2955.0, 3012.0, 3078.0, 3142.0, 3226.0, 3428.0, 3634.0, 3802.0, 3955.0, 4078.0, 4192.0, 4290.0, 4377.0, 4456.0, 4529.0, 4596.0, 4658.0]),
                       list([2236.0, 2452.0, 2571.0, 2707.0, 2784.0, 2829.0, 2867.0, 2902.0, 2907.0, 2914.0, 2924.0, 2934.0, 2972.0, 2996.0, 3034.0, 3076.0, 3126.0, 3202.0, 3428.0, 3606.0, 3763.0, 3906.0, 4033.0, 4153.0, 4259.0, 4350.0, 4432.0, 4508.0, 4581.0, 4653.0]),
                       list([2142.0, 2406.0, 2545.0, 2706.0, 2784.0, 2831.0, 2870.0, 2906.0, 2911.0, 2919.0, 2931.0, 2956.0, 2998.0, 3016.0, 3052.0, 3088.0, 3127.0, 3196.0, 3404.0, 3582.0, 3731.0, 3868.0, 3998.0, 4120.0, 4232.0, 4330.0, 4421.0, 4505.0, 4588.0, 4671.0]),
                       list([2053.0, 2362.0, 2535.0, 2705.0, 2785.0, 2834.0, 2875.0, 2914.0, 2921.0, 2931.0, 2946.0, 2976.0, 3013.0, 3035.0, 3075.0, 3113.0, 3145.0, 3206.0, 3389.0, 3560.0, 3704.0, 3842.0, 3976.0, 4096.0, 4212.0, 4315.0, 4417.0, 4513.0, 4609.0, 4710.0]),
                       list([1961.0, 2330.0, 2531.0, 2705.0, 2785.0, 2841.0, 2890.0, 2928.0, 2934.0, 2944.0, 2959.0, 2992.0, 3027.0, 3052.0, 3098.0, 3136.0, 3169.0, 3220.0, 3379.0, 3546.0, 3696.0, 3830.0, 3960.0, 4081.0, 4198.0, 4310.0, 4422.0, 4531.0, 4644.0, 4771.0]),
                       list([1714.0, 2280.0, 2525.0, 2699.0, 2784.0, 2849.0, 2903.0, 2948.0, 2956.0, 2968.0, 2986.0, 3017.0, 3045.0, 3069.0, 3115.0, 3155.0, 3192.0, 3234.0, 3364.0, 3518.0, 3669.0, 3798.0, 3923.0, 4049.0, 4184.0, 4327.0, 4481.0, 4647.0, 4837.0, 5064.0]),
                       list([1570.0, 2261.0, 2519.0, 2682.0, 2768.0, 2843.0, 2900.0, 2946.0, 2954.0, 2967.0, 2987.0, 3021.0, 3050.0, 3075.0, 3118.0, 3159.0, 3197.0, 3240.0, 3363.0, 3508.0, 3651.0, 3781.0, 3908.0, 4041.0, 4191.0, 4363.0, 4553.0, 4768.0, 5038.0, 5363.0]),
                       list([2237.0, 2495.0, 2614.0, 2721.0, 2827.0, 2892.0, 2943.0, 2952.0, 2965.0, 2986.0, 3023.0, 3056.0, 3086.0, 3131.0, 3169.0, 3205.0, 3252.0, 3374.0, 3498.0, 3621.0, 3748.0, 3891.0, 4054.0, 4263.0, 4530.0, 4850.0, 5264.0, 5754.0, 6422.0]),
                       list([2209.0, 2411.0, 2555.0, 2681.0, 2811.0, 2883.0, 2940.0, 2950.0, 2964.0, 2986.0, 3025.0, 3059.0, 3090.0, 3141.0, 3182.0, 3215.0, 3263.0, 3392.0, 3500.0, 3609.0, 3733.0, 3894.0, 4102.0, 4406.0, 4793.0, 5302.0, 5911.0, 6330.0, 6475.0]),
                       list([2181.0, 2306.0, 2500.0, 2655.0, 2792.0, 2874.0, 2935.0, 2946.0, 2961.0, 2985.0, 3026.0, 3061.0, 3094.0, 3148.0, 3194.0, 3228.0, 3276.0, 3410.0, 3507.0, 3604.0, 3731.0, 3916.0, 4186.0, 4601.0, 5147.0, 5796.0, 6102.0, 6271.0, 6622.0]),
                       list([2139.0, 2184.0, 2448.0, 2632.0, 2755.0, 2860.0, 2925.0, 2936.0, 2952.0, 2977.0, 3020.0, 3056.0, 3089.0, 3145.0, 3193.0, 3231.0, 3279.0, 3413.0, 3508.0, 3602.0, 3739.0, 3960.0, 4320.0, 4848.0, 5520.0, 5897.0, 6052.0, 6335.0, 6743.0]),
                       list([1996.0, 2005.0, 2366.0, 2581.0, 2699.0, 2830.0, 2904.0, 2916.0, 2933.0, 2961.0, 3007.0, 3046.0, 3081.0, 3138.0, 3187.0, 3229.0, 3280.0, 3410.0, 3509.0, 3606.0, 3771.0, 4084.0, 4645.0, 5301.0, 5570.0, 5822.0, 6090.0, 6413.0, 6757.0]),
                       list([1760.0, 1838.0, 2304.0, 2525.0, 2671.0, 2804.0, 2882.0, 2896.0, 2916.0, 2945.0, 2996.0, 3037.0, 3073.0, 3133.0, 3182.0, 3225.0, 3280.0, 3407.0, 3512.0, 3617.0, 3821.0, 4270.0, 4938.0, 5236.0, 5536.0, 5833.0, 6133.0, 6435.0, 6759.0]),
                       list([2047.0, 2387.0, 2586.0, 2709.0, 2829.0, 2845.0, 2869.0, 2905.0, 2962.0, 3010.0, 3052.0, 3120.0, 3173.0, 3217.0, 3275.0, 3406.0, 3522.0, 3668.0, 4039.0, 4497.0, 4830.0, 5212.0, 5533.0, 5876.0, 6164.0, 6451.0, 6766.0]),
                       list([1912.0, 2311.0, 2525.0, 2670.0, 2794.0, 2813.0, 2839.0, 2879.0, 2942.0, 2993.0, 3038.0, 3112.0, 3167.0, 3214.0, 3272.0, 3409.0, 3529.0, 3705.0, 4100.0, 4402.0, 4826.0, 5209.0, 5541.0, 5886.0, 6170.0, 6453.0, 6766.0]),
                       list([1779.0, 2225.0, 2471.0, 2635.0, 2757.0, 2780.0, 2810.0, 2853.0, 2924.0, 2978.0, 3026.0, 3105.0, 3162.0, 3210.0, 3270.0, 3411.0, 3536.0, 3727.0, 4030.0, 4400.0, 4825.0, 5208.0, 5547.0, 5889.0, 6172.0, 6455.0, 6767.0]),
                       list([1940.0, 2268.0, 2488.0, 2639.0, 2664.0, 2698.0, 2756.0, 2853.0, 2926.0, 2982.0, 3076.0, 3144.0, 3198.0, 3264.0, 3414.0, 3537.0, 3701.0, 3981.0, 4390.0, 4824.0, 5207.0, 5567.0, 5895.0, 6177.0, 6456.0, 6766.0]),
                       list([1667.0, 2063.0, 2331.0, 2519.0, 2552.0, 2594.0, 2660.0, 2782.0, 2874.0, 2946.0, 3055.0, 3134.0, 3192.0, 3261.0, 3415.0, 3528.0, 3689.0, 3971.0, 4388.0, 4822.0, 5209.0, 5572.0, 5898.0, 6179.0, 6456.0, 6761.0]),
                       list([1836.0, 2207.0, 2419.0, 2457.0, 2508.0, 2594.0, 2726.0, 2844.0, 2929.0, 3048.0, 3131.0, 3191.0, 3260.0, 3415.0, 3526.0, 3685.0, 3969.0, 4387.0, 4821.0, 5210.0, 5574.0, 5901.0, 6181.0, 6453.0, 6750.0]),
                       list([1698.0, 2082.0, 2335.0, 2378.0, 2435.0, 2539.0, 2696.0, 2828.0, 2921.0, 3046.0, 3130.0, 3191.0, 3260.0, 3415.0, 3524.0, 3682.0, 3968.0, 4386.0, 4821.0, 5211.0, 5576.0, 5904.0, 6182.0, 6451.0, 6742.0]),
                       list([1555.0, 1908.0, 2248.0, 2301.0, 2370.0, 2485.0, 2671.0, 2817.0, 2917.0, 3045.0, 3130.0, 3191.0, 3260.0, 3415.0, 3522.0, 3681.0, 3969.0, 4385.0, 4821.0, 5213.0, 5580.0, 5907.0, 6183.0, 6446.0, 6727.0]),
                       list([1716.0, 2136.0, 2212.0, 2294.0, 2436.0, 2656.0, 2812.0, 2915.0, 3045.0, 3130.0, 3191.0, 3260.0, 3415.0, 3521.0, 3680.0, 3968.0, 4386.0, 4823.0, 5216.0, 5586.0, 5910.0, 6184.0, 6442.0, 6705.0]),
                       list([1584.0, 2047.0, 2109.0, 2243.0, 2400.0, 2648.0, 2810.0, 2915.0, 3045.0, 3130.0, 3191.0, 3261.0, 3415.0, 3520.0, 3680.0, 3967.0, 4387.0, 4826.0, 5220.0, 5592.0, 5916.0, 6187.0, 6438.0, 6684.0]),
                       list([1720.0, 1811.0, 2083.0, 2348.0, 2643.0, 2810.0, 2915.0, 3045.0, 3130.0, 3191.0, 3261.0, 3415.0, 3519.0, 3679.0, 3969.0, 4396.0, 4842.0, 5244.0, 5627.0, 5945.0, 6202.0, 6403.0, 6464.0]),
                       list([1650.0, 1793.0, 2052.0, 2343.0, 2644.0, 2810.0, 2915.0, 3046.0, 3131.0, 3191.0, 3261.0, 3416.0, 3520.0, 3679.0, 3974.0, 4407.0, 4859.0, 5270.0, 5661.0, 5974.0, 6191.0, 6177.0]),
                       list([1631.0, 1790.0, 2041.0, 2343.0, 2644.0, 2811.0, 2916.0, 3046.0, 3131.0, 3192.0, 3262.0, 3416.0, 3520.0, 3680.0, 3979.0, 4418.0, 4878.0, 5297.0, 5697.0, 5996.0, 6053.0]),
                       list([1627.0, 1790.0, 2038.0, 2343.0, 2644.0, 2811.0, 2916.0, 3046.0, 3131.0, 3192.0, 3262.0, 3416.0, 3521.0, 3682.0, 3984.0, 4430.0, 4896.0, 5325.0, 5732.0, 5971.0]),
                       list([1626.0, 1789.0, 2036.0, 2344.0, 2645.0, 2812.0, 2917.0, 3047.0, 3132.0, 3193.0, 3264.0, 3416.0, 3522.0, 3686.0, 4001.0, 4467.0, 4958.0, 5420.0, 5772.0]),
                       list([1626.0, 1790.0, 2037.0, 2345.0, 2646.0, 2812.0, 2917.0, 3048.0, 3133.0, 3194.0, 3265.0, 3416.0, 3524.0, 3689.0, 4013.0, 4493.0, 5002.0, 5495.0])], dtype='object')
    radii = np.array([list([235842300.00000003, 287324100.0, 324196200.0, 416724300.0, 459857700.0, 550298700.0, 585083700.0, 626130000.0, 651870900.0, 672046200.0, 723528000.0, 765965700.0, 718658100.0, 775705500.0, 883539000.0, 982328400.0, 1088770500.0, 1202169600.0, 1532627100.0, 1626546600.0, 2337552000.0, 2372337000.0, 2386251000.0, 2393208000.0, 2400165000.0, 2407122000.0, 2414079000.0, 2414079000.0, 2414079000.0, 2414079000.0]),
                      list([188534700.0, 226798200.0, 258800400.0, 324891900.0, 372895200.0, 436899600.0, 488381400.0, 543341700.0, 558647100.0, 576735300.0, 610128900.0, 667176300.0, 697091400.0, 750660300.0, 844579800.0, 923193900.0, 1019200500.0, 1045637099.9999999, 1138165200.0, 1219562100.0, 1300959000.0, 1371224700.0, 2532348000.0, 2546262000.0, 2546262000.0, 2553219000.0, 2560176000.0, 2560176000.0, 2567133000.0, 2567133000.0]),
                      list([153054000.0, 182273400.0, 207318600.0, 269235900.0, 334631700.0, 399331799.99999994, 459857700.0, 516905100.0, 529427700.0, 546124500.0, 575343900.0, 628912800.0, 619173000.0, 641435400.0, 667176300.0, 688743000.0, 737442000.0, 749268900.0, 861276600.0, 948239100.0, 1026853200.0, 1086683400.0, 1154166300.0, 1216083600.0, 1275218100.0, 2692359000.0, 2699316000.0, 2699316000.0, 2706273000.0, 2706273000.0]),
                      list([137748600.0, 161402400.0, 185751900.0, 259496100.0, 326979000.0, 388896300.00000006, 441073800.0, 483511499.99999994, 493947000.0, 505773900.0, 527340600.0, 552385800.0, 493947000.0, 512035200.0, 544733100.0, 575343900.0, 613607400.0, 637956900.0, 742311900.0, 825100200.0, 898844400.0, 955891800.0000001, 1014330600.0, 1067899500.0, 1115902800.0, 1164601800.0, 2775843000.0, 2782800000.0, 2782800000.0, 2782800000.0]),
                      list([128008800.0, 149575500.0, 178794900.0, 255321900.0, 317934900.0, 371503800.0, 404897400.0, 413941500.0, 418115700.0, 421594200.0, 434812500.0, 450117900.0, 427159800.0, 444552300.0, 478641599.99999994, 507165300.0, 544037400.0, 571865400.0, 674133300.0, 749964600.0, 816056100.0, 875190600.0, 928063800.0, 973979999.9999999, 1019200500.0, 1062333899.9999999, 1113815700.0, 2831499000.0, 2831499000.0, 2831499000.0]),
                      list([122443200.0, 141922800.0, 176707800.0, 251843400.0, 306108000.0, 340197300.0, 345067200.0, 355502700.0, 360372600.0, 365242500.0, 376373700.0, 393766199.99999994, 385417800.00000006, 403506000.0, 436203900.0, 465423300.0, 496034100.0, 527340600.0, 625434300.0, 699874200.0, 761791500.0, 816056100.0, 864755100.0000001, 908584200.0, 950326200.0000001, 994155300.0, 1039375800.0, 1091553300.0, 2866284000.0, 2859327000.0]),
                      list([112007700.0, 133574400.0, 172533600.0, 228885300.0, 231668100.0, 245582100.0, 260887500.0, 275497200.0, 278975700.0, 285932700.0, 294976800.0, 309586500.0, 315152100.0, 331848900.0, 359676900.0, 385417800.00000006, 412550100.0, 443160900.0, 532210500.0, 601084800.0, 656740800.0, 702657000.0, 745094700.0, 788228100.0, 829274400.0, 873799200.0, 919019700.0, 972588599.9999999, 2928897000.0, 2901069000.0]),
                      list([108529200.0, 131487300.0, 169750800.0, 196883099.99999997, 203144400.0, 219145500.0, 233059500.0, 249060600.0, 252539100.0, 256713300.0, 266453100.0, 281758500.0, 289411200.0, 302629500.0, 330457500.0, 354111300.0, 378460800.0, 408375900.0, 496034100.0, 560038500.0, 612216000.0, 655349400.0, 699178499.9999999, 740224800.0, 779879700.0, 827187300.0, 875886299.9999999, 937107900.0, 2935854000.0, 2894112000.0]),
                      list([128704500.0, 159315300.0, 154445400.0, 166968000.0, 181577700.0, 196187399.99999997, 210101400.0, 213579900.0, 217058400.0, 225406800.0, 237929400.00000003, 246973500.0, 259496100.0, 283149900.0, 304716600.0, 326283300.0, 353415600.0, 433421100.0, 491859900.0, 538471800.0, 580213800.0, 622651500.0, 664393500.0, 709614000.0, 759704400.0, 827187300.0, 2970638999.9999995, 2887155000.0000005, 2831499000.0]),
                      list([125921700.0, 137748600.0, 134965800.0, 148184100.0, 162098100.0, 174620700.0, 188534700.0, 191317500.00000003, 194100300.00000003, 201753000.0, 212884200.0, 221928300.0, 233755200.0, 254626200.0, 275497200.0, 293585400.0, 320022000.0, 393070499.99999994, 447335100.0, 493947000.0, 536384700.0, 576735300.0, 622651500.0, 672046200.0, 735354900.0, 3005424000.0, 2928897000.0, 2991510000.0, 2977596000.0]),
                      list([123138900.0, 119660399.99999999, 123834600.0, 136357200.0, 148879800.0, 160706700.0, 173229300.0, 176012100.0, 178794900.0, 185056200.0, 196883099.99999997, 205231500.0, 215667000.0, 235146600.00000003, 253930500.0, 272018700.0, 295672500.0, 362459700.0, 416724300.0, 460553400.0, 504382500.0, 548907300.0, 596910600.0, 652566600.0, 734659200.0, 2998466999.9999995, 3054123000.0, 3019338000.0, 2970638999.9999995]),
                      list([118964700.00000001, 110616300.0, 116181900.0, 128704500.0, 139140000.0, 150966900.0, 162098100.0, 164185200.0, 167663700.0, 173229300.0, 183664800.0, 192708900.00000003, 201057300.0, 220536900.0, 237929400.00000003, 254626200.0, 276888600.0, 339501600.0, 392374799.99999994, 436203900.0, 481424399.99999994, 527340600.0, 580213800.0, 643522500.0, 732572100.0, 3088908000.0000005, 3068037000.0, 3012381000.0, 2984553000.0]),
                      list([109920600.0, 100180799.99999999, 106442100.0, 116181900.0, 127313100.0, 136357200.0, 146792700.0, 148879800.0, 151662600.0, 156532500.0, 166272300.0, 173925000.0, 182969100.0, 199665899.99999997, 214971300.0, 229581000.0, 250452000.0, 308195100.0, 357589800.0, 403506000.0, 450117900.0, 502991100.0, 560734200.0, 620564400.0, 626130000.0, 3109779000.0, 3061080000.0000005, 3012381000.0, 2991510000.0]),
                      list([101572200.0, 95310900.0, 100180799.99999999, 108529200.0, 118269000.00000001, 126617400.0, 135661500.0, 137748600.0, 140531400.0, 145401300.0, 153749700.0, 161402400.0, 169750800.0, 185056200.0, 198274499.99999997, 212188500.0, 230972400.0, 285932700.0, 333240300.0, 381243600.0, 432029700.0, 488381400.0, 539863200.0, 557951400.0, 622651500.0, 3109779000.0, 3061080000.0000005, 3019338000.0, 2991510000.0]),
                      list([89745300.0, 96006600.00000001, 102963600.0, 109920600.0, 117573300.00000001, 118964700.00000001, 121051799.99999999, 124530300.0, 132183000.0, 137748600.0, 144705600.0, 157923900.0, 170446500.0, 181577700.0, 198274499.99999997, 246973500.0, 292889700.0, 345067200.0, 405593100.0, 442465200.0, 490468500.0, 555864300.0, 626825700.0, 3102822000.0, 3054123000.0, 3012381000.0, 2991510000.0]),
                      list([86266800.0, 90441000.0, 96702300.00000001, 102963600.0, 109920600.0, 111312000.0, 112703400.0, 116181900.0, 122443200.0, 128704500.0, 134965800.0, 146792700.0, 157923900.0, 169055100.0, 184360500.0, 231668100.0, 277584300.0, 333936000.0, 392374799.99999994, 434116800.0, 491859900.0, 557255700.0, 628217100.0, 3102822000.0, 3054123000.0, 3012381000.0, 2991510000.0]),
                      list([83484000.0, 86962500.0, 92528100.0, 98093699.99999999, 103659300.0, 105050700.0, 107137800.0, 109920600.0, 116181900.0, 121747499.99999999, 127313100.0, 138444300.0, 149575500.0, 159315300.0, 174620700.0, 221232600.0, 267148800.0, 327674700.0, 379852200.0, 435508200.0, 492555600.0, 557255700.0, 628912800.0, 3102822000.0, 3054123000.0, 3012381000.0, 2984553000.0]),
                      list([79309800.0, 82092600.0, 85571100.0, 89745300.0, 91136700.0, 91832400.0, 94615200.0, 99485099.99999999, 105050700.0, 109920600.0, 120356099.99999999, 130791600.0, 141227100.0, 157228200.0, 204535800.0, 253930500.0, 315152100.0, 380547900.0, 438291000.0, 494642700.0, 559342800.0, 630999900.0, 3102822000.0, 3054123000.0, 3012381000.0, 2984553000.0]),
                      list([75135600.0, 75831300.0, 77918400.0, 81396900.0, 82092600.0, 82788300.0, 84875400.0, 89745300.0, 95310900.0, 100180799.99999999, 112007700.0, 123834600.0, 134965800.0, 151662600.0, 201753000.0, 253234800.0, 316543500.0, 382635000.00000006, 439682400.0, 496034100.0, 560734200.0, 632391300.0, 3102822000.0, 3054123000.0, 3005424000.0, 2977596000.0]),
                      list([72352800.0, 73744200.0, 76527000.0, 77222700.0, 77918400.0, 80005500.0, 84875400.0, 90441000.0, 96702300.00000001, 109224900.0, 121747499.99999999, 133574400.0, 150966900.0, 201753000.0, 254626200.0, 317934900.0, 383330700.00000006, 440378100.0, 496729800.0, 562125600.0, 633782700.0, 3095865000.0, 3047166000.0, 3005424000.0, 2970638999.9999995]),
                      list([70265700.0, 70961400.0, 73048500.0, 73744200.0, 74439900.0, 77222700.0, 82092600.0, 88353900.0, 95310900.0, 108529200.0, 121747499.99999999, 133574400.0, 150966900.0, 201753000.0, 255321900.0, 317934900.0, 384026400.00000006, 441073800.0, 497425500.0, 563517000.0, 635174100.0, 3095865000.0, 3047166000.0, 2998466999.9999995, 2963682000.0]),
                      list([68874300.0, 68178600.0, 70265700.0, 70961400.0, 71657100.0, 74439900.0, 80005500.0, 87658200.0, 94615200.0, 108529200.0, 121051799.99999999, 133574400.0, 150966900.0, 201753000.0, 256713300.0, 318630600.0, 384026400.00000006, 441769500.0, 498121200.0, 564212700.0, 637261200.0, 3088908000.0000005, 3040209000.0, 2998466999.9999995, 2956725000.0]),
                      list([66787200.0, 67482900.0, 68178600.0, 69570000.0, 72352800.0, 79309800.0, 86962500.0, 93919500.0, 108529200.0, 121051799.99999999, 133574400.0, 150966900.0, 201753000.0, 258104700.0, 319326300.0, 384722100.00000006, 442465200.0, 499512600.0, 565604100.0, 640044000.0, 3088908000.0000005, 3033252000.0, 2984553000.0, 2949768000.0]),
                      list([65395800.0, 65395800.0, 66091500.0, 67482900.0, 70961400.0, 78614100.0, 86266800.0, 93919500.0, 108529200.0, 121051799.99999999, 133574400.0, 150966900.0, 202448700.0, 258800400.0, 320022000.0, 385417800.00000006, 443160900.0, 500208300.0, 567691200.0, 642131100.0, 3081951000.0, 3033252000.0, 2977596000.0, 2935854000.0]),
                      list([63308700.0, 63308700.0, 64700100.0, 68874300.0, 78614100.0, 86266800.0, 93919500.0, 108529200.0, 121051799.99999999, 133574400.0, 150966900.0, 203840100.0, 260191800.0, 321413400.0, 387504900.00000006, 445248000.0, 504382500.0, 575343900.0, 656740800.0, 3061080000.0000005, 2998466999.9999995, 2928897000.0, 2852369999.9999995]),
                      list([62613000.0, 63308700.0, 64004400.0, 68874300.0, 78614100.0, 86266800.0, 93919500.0, 108529200.0, 121051799.99999999, 133574400.0, 150966900.0, 205231500.0, 260887500.0, 322804800.0, 388896300.00000006, 447335100.0, 508556700.0, 582996600.0, 672046200.0, 3033252000.0, 2956725000.0, 2859327000.0]),
                      list([62613000.0, 63308700.0, 64004400.0, 68874300.0, 78614100.0, 86266800.0, 93919500.0, 108529200.0, 121051799.99999999, 133574400.0, 150966900.0, 206622900.0, 261583200.0, 323500500.0, 390287700.00000006, 449422200.0, 512035200.0, 591345000.0, 691525800.0, 3005424000.0, 2901069000.0]),
                      list([62613000.0, 63308700.0, 63308700.0, 68874300.0, 78614100.0, 86266800.0, 93919500.0, 108529200.0, 121747499.99999999, 133574400.0, 151662600.0, 208014300.0, 262278900.0, 324891900.0, 390983400.00000006, 450813600.0, 516209400.0, 600389100.0, 3068037000.0, 2963682000.0]),
                      list([62613000.0, 63308700.0, 63308700.0, 68874300.0, 78614100.0, 86266800.0, 93919500.0, 108529200.0, 121747499.99999999, 134270100.0, 152358300.0, 210797100.0, 263670300.0, 326283300.0, 393766199.99999994, 456379200.0, 529427700.0, 635869800.0, 2984553000.0]),
                      list([62613000.0, 63308700.0, 63308700.0, 68874300.0, 78614100.0, 86266800.0, 93919500.0, 108529200.0, 121747499.99999999, 134270100.0, 153054000.0, 211492800.0, 264366000.0, 327674700.0, 395853299.99999994, 460553400.0, 539863200.0, 670654800.0])], dtype='object')
    
    masses = np.array([list([1.988435e+28, 2.9826525e+28, 3.97687e+28, 5.965305e+28, 7.95374e+28, 9.942175000000001e+28, 1.193061e+29, 1.3919045000000002e+29, 1.4316731999999999e+29, 1.4913262499999998e+29, 1.590748e+29, 1.7895915e+29, 1.9884350000000003e+29, 2.1872785e+29, 2.5849655e+29, 2.9826524999999997e+29, 3.3803395000000004e+29, 3.9768700000000005e+29, 5.9653049999999994e+29, 7.953740000000001e+29, 9.942175e+29, 1.1930609999999999e+30, 1.3919044999999998e+30, 1.5907480000000002e+30, 1.7895915e+30, 1.988435e+30, 2.1872785000000001e+30, 2.3861219999999998e+30, 2.5849655e+30, 2.7838089999999995e+30]),
                        list([1.988435e+28, 2.9826525e+28, 3.97687e+28, 5.965305e+28, 7.95374e+28, 9.942175000000001e+28, 1.193061e+29, 1.3919045000000002e+29, 1.4316731999999999e+29, 1.4913262499999998e+29, 1.590748e+29, 1.7895915e+29, 1.9884350000000003e+29, 2.1872785e+29, 2.5849655e+29, 2.9826524999999997e+29, 3.3803395000000004e+29, 3.9768700000000005e+29, 5.9653049999999994e+29, 7.953740000000001e+29, 9.942175e+29, 1.1930609999999999e+30, 1.3919044999999998e+30, 1.5907480000000002e+30, 1.7895915e+30, 1.988435e+30, 2.1872785000000001e+30, 2.3861219999999998e+30, 2.5849655e+30, 2.7838089999999995e+30]),
                        list([1.988435e+28, 2.9826525e+28, 3.97687e+28, 5.965305e+28, 7.95374e+28, 9.942175000000001e+28, 1.193061e+29, 1.3919045000000002e+29, 1.4316731999999999e+29, 1.4913262499999998e+29, 1.590748e+29, 1.7895915e+29, 1.9884350000000003e+29, 2.1872785e+29, 2.5849655e+29, 2.9826524999999997e+29, 3.3803395000000004e+29, 3.9768700000000005e+29, 5.9653049999999994e+29, 7.953740000000001e+29, 9.942175e+29, 1.1930609999999999e+30, 1.3919044999999998e+30, 1.5907480000000002e+30, 1.7895915e+30, 1.988435e+30, 2.1872785000000001e+30, 2.3861219999999998e+30, 2.5849655e+30, 2.7838089999999995e+30]),
                        list([1.988435e+28, 2.9826525e+28, 3.97687e+28, 5.965305e+28, 7.95374e+28, 9.942175000000001e+28, 1.193061e+29, 1.3919045000000002e+29, 1.4316731999999999e+29, 1.4913262499999998e+29, 1.590748e+29, 1.7895915e+29, 1.9884350000000003e+29, 2.1872785e+29, 2.5849655e+29, 2.9826524999999997e+29, 3.3803395000000004e+29, 3.9768700000000005e+29, 5.9653049999999994e+29, 7.953740000000001e+29, 9.942175e+29, 1.1930609999999999e+30, 1.3919044999999998e+30, 1.5907480000000002e+30, 1.7895915e+30, 1.988435e+30, 2.1872785000000001e+30, 2.3861219999999998e+30, 2.5849655e+30, 2.7838089999999995e+30]),
                        list([1.988435e+28, 2.9826525e+28, 3.97687e+28, 5.965305e+28, 7.95374e+28, 9.942175000000001e+28, 1.193061e+29, 1.3919045000000002e+29, 1.4316731999999999e+29, 1.4913262499999998e+29, 1.590748e+29, 1.7895915e+29, 1.9884350000000003e+29, 2.1872785e+29, 2.5849655e+29, 2.9826524999999997e+29, 3.3803395000000004e+29, 3.9768700000000005e+29, 5.9653049999999994e+29, 7.953740000000001e+29, 9.942175e+29, 1.1930609999999999e+30, 1.3919044999999998e+30, 1.5907480000000002e+30, 1.7895915e+30, 1.988435e+30, 2.1872785000000001e+30, 2.3861219999999998e+30, 2.5849655e+30, 2.7838089999999995e+30]),
                        list([1.988435e+28, 2.9826525e+28, 3.97687e+28, 5.965305e+28, 7.95374e+28, 9.942175000000001e+28, 1.193061e+29, 1.3919045000000002e+29, 1.4316731999999999e+29, 1.4913262499999998e+29, 1.590748e+29, 1.7895915e+29, 1.9884350000000003e+29, 2.1872785e+29, 2.5849655e+29, 2.9826524999999997e+29, 3.3803395000000004e+29, 3.9768700000000005e+29, 5.9653049999999994e+29, 7.953740000000001e+29, 9.942175e+29, 1.1930609999999999e+30, 1.3919044999999998e+30, 1.5907480000000002e+30, 1.7895915e+30, 1.988435e+30, 2.1872785000000001e+30, 2.3861219999999998e+30, 2.5849655e+30, 2.7838089999999995e+30]),
                        list([1.988435e+28, 2.9826525e+28, 3.97687e+28, 5.965305e+28, 7.95374e+28, 9.942175000000001e+28, 1.193061e+29, 1.3919045000000002e+29, 1.4316731999999999e+29, 1.4913262499999998e+29, 1.590748e+29, 1.7895915e+29, 1.9884350000000003e+29, 2.1872785e+29, 2.5849655e+29, 2.9826524999999997e+29, 3.3803395000000004e+29, 3.9768700000000005e+29, 5.9653049999999994e+29, 7.953740000000001e+29, 9.942175e+29, 1.1930609999999999e+30, 1.3919044999999998e+30, 1.5907480000000002e+30, 1.7895915e+30, 1.988435e+30, 2.1872785000000001e+30, 2.3861219999999998e+30, 2.5849655e+30, 2.7838089999999995e+30]),
                        list([1.988435e+28, 2.9826525e+28, 3.97687e+28, 5.965305e+28, 7.95374e+28, 9.942175000000001e+28, 1.193061e+29, 1.3919045000000002e+29, 1.4316731999999999e+29, 1.4913262499999998e+29, 1.590748e+29, 1.7895915e+29, 1.9884350000000003e+29, 2.1872785e+29, 2.5849655e+29, 2.9826524999999997e+29, 3.3803395000000004e+29, 3.9768700000000005e+29, 5.9653049999999994e+29, 7.953740000000001e+29, 9.942175e+29, 1.1930609999999999e+30, 1.3919044999999998e+30, 1.5907480000000002e+30, 1.7895915e+30, 1.988435e+30, 2.1872785000000001e+30, 2.3861219999999998e+30, 2.5849655e+30, 2.7838089999999995e+30]),
                        list([2.9826525e+28, 3.97687e+28, 5.965305e+28, 7.95374e+28, 9.942175000000001e+28, 1.193061e+29, 1.3919045000000002e+29, 1.4316731999999999e+29, 1.4913262499999998e+29, 1.590748e+29, 1.7895915e+29, 1.9884350000000003e+29, 2.1872785e+29, 2.5849655e+29, 2.9826524999999997e+29, 3.3803395000000004e+29, 3.9768700000000005e+29, 5.9653049999999994e+29, 7.953740000000001e+29, 9.942175e+29, 1.1930609999999999e+30, 1.3919044999999998e+30, 1.5907480000000002e+30, 1.7895915e+30, 1.988435e+30, 2.1872785000000001e+30, 2.3861219999999998e+30, 2.5849655e+30, 2.7838089999999995e+30]),
                        list([2.9826525e+28, 3.97687e+28, 5.965305e+28, 7.95374e+28, 9.942175000000001e+28, 1.193061e+29, 1.3919045000000002e+29, 1.4316731999999999e+29, 1.4913262499999998e+29, 1.590748e+29, 1.7895915e+29, 1.9884350000000003e+29, 2.1872785e+29, 2.5849655e+29, 2.9826524999999997e+29, 3.3803395000000004e+29, 3.9768700000000005e+29, 5.9653049999999994e+29, 7.953740000000001e+29, 9.942175e+29, 1.1930609999999999e+30, 1.3919044999999998e+30, 1.5907480000000002e+30, 1.7895915e+30, 1.988435e+30, 2.1872785000000001e+30, 2.3861219999999998e+30, 2.5849655e+30, 2.7838089999999995e+30]),
                        list([2.9826525e+28, 3.97687e+28, 5.965305e+28, 7.95374e+28, 9.942175000000001e+28, 1.193061e+29, 1.3919045000000002e+29, 1.4316731999999999e+29, 1.4913262499999998e+29, 1.590748e+29, 1.7895915e+29, 1.9884350000000003e+29, 2.1872785e+29, 2.5849655e+29, 2.9826524999999997e+29, 3.3803395000000004e+29, 3.9768700000000005e+29, 5.9653049999999994e+29, 7.953740000000001e+29, 9.942175e+29, 1.1930609999999999e+30, 1.3919044999999998e+30, 1.5907480000000002e+30, 1.7895915e+30, 1.988435e+30, 2.1872785000000001e+30, 2.3861219999999998e+30, 2.5849655e+30, 2.7838089999999995e+30]),
                        list([2.9826525e+28, 3.97687e+28, 5.965305e+28, 7.95374e+28, 9.942175000000001e+28, 1.193061e+29, 1.3919045000000002e+29, 1.4316731999999999e+29, 1.4913262499999998e+29, 1.590748e+29, 1.7895915e+29, 1.9884350000000003e+29, 2.1872785e+29, 2.5849655e+29, 2.9826524999999997e+29, 3.3803395000000004e+29, 3.9768700000000005e+29, 5.9653049999999994e+29, 7.953740000000001e+29, 9.942175e+29, 1.1930609999999999e+30, 1.3919044999999998e+30, 1.5907480000000002e+30, 1.7895915e+30, 1.988435e+30, 2.1872785000000001e+30, 2.3861219999999998e+30, 2.5849655e+30, 2.7838089999999995e+30]),
                        list([2.9826525e+28, 3.97687e+28, 5.965305e+28, 7.95374e+28, 9.942175000000001e+28, 1.193061e+29, 1.3919045000000002e+29, 1.4316731999999999e+29, 1.4913262499999998e+29, 1.590748e+29, 1.7895915e+29, 1.9884350000000003e+29, 2.1872785e+29, 2.5849655e+29, 2.9826524999999997e+29, 3.3803395000000004e+29, 3.9768700000000005e+29, 5.9653049999999994e+29, 7.953740000000001e+29, 9.942175e+29, 1.1930609999999999e+30, 1.3919044999999998e+30, 1.5907480000000002e+30, 1.7895915e+30, 1.988435e+30, 2.1872785000000001e+30, 2.3861219999999998e+30, 2.5849655e+30, 2.7838089999999995e+30]),
                        list([2.9826525e+28, 3.97687e+28, 5.965305e+28, 7.95374e+28, 9.942175000000001e+28, 1.193061e+29, 1.3919045000000002e+29, 1.4316731999999999e+29, 1.4913262499999998e+29, 1.590748e+29, 1.7895915e+29, 1.9884350000000003e+29, 2.1872785e+29, 2.5849655e+29, 2.9826524999999997e+29, 3.3803395000000004e+29, 3.9768700000000005e+29, 5.9653049999999994e+29, 7.953740000000001e+29, 9.942175e+29, 1.1930609999999999e+30, 1.3919044999999998e+30, 1.5907480000000002e+30, 1.7895915e+30, 1.988435e+30, 2.1872785000000001e+30, 2.3861219999999998e+30, 2.5849655e+30, 2.7838089999999995e+30]),
                        list([5.965305e+28, 7.95374e+28, 9.942175000000001e+28, 1.193061e+29, 1.3919045000000002e+29, 1.4316731999999999e+29, 1.4913262499999998e+29, 1.590748e+29, 1.7895915e+29, 1.9884350000000003e+29, 2.1872785e+29, 2.5849655e+29, 2.9826524999999997e+29, 3.3803395000000004e+29, 3.9768700000000005e+29, 5.9653049999999994e+29, 7.953740000000001e+29, 9.942175e+29, 1.1930609999999999e+30, 1.3919044999999998e+30, 1.5907480000000002e+30, 1.7895915e+30, 1.988435e+30, 2.1872785000000001e+30, 2.3861219999999998e+30, 2.5849655e+30, 2.7838089999999995e+30]),
                        list([5.965305e+28, 7.95374e+28, 9.942175000000001e+28, 1.193061e+29, 1.3919045000000002e+29, 1.4316731999999999e+29, 1.4913262499999998e+29, 1.590748e+29, 1.7895915e+29, 1.9884350000000003e+29, 2.1872785e+29, 2.5849655e+29, 2.9826524999999997e+29, 3.3803395000000004e+29, 3.9768700000000005e+29, 5.9653049999999994e+29, 7.953740000000001e+29, 9.942175e+29, 1.1930609999999999e+30, 1.3919044999999998e+30, 1.5907480000000002e+30, 1.7895915e+30, 1.988435e+30, 2.1872785000000001e+30, 2.3861219999999998e+30, 2.5849655e+30, 2.7838089999999995e+30]),
                        list([5.965305e+28, 7.95374e+28, 9.942175000000001e+28, 1.193061e+29, 1.3919045000000002e+29, 1.4316731999999999e+29, 1.4913262499999998e+29, 1.590748e+29, 1.7895915e+29, 1.9884350000000003e+29, 2.1872785e+29, 2.5849655e+29, 2.9826524999999997e+29, 3.3803395000000004e+29, 3.9768700000000005e+29, 5.9653049999999994e+29, 7.953740000000001e+29, 9.942175e+29, 1.1930609999999999e+30, 1.3919044999999998e+30, 1.5907480000000002e+30, 1.7895915e+30, 1.988435e+30, 2.1872785000000001e+30, 2.3861219999999998e+30, 2.5849655e+30, 2.7838089999999995e+30]),
                        list([7.95374e+28, 9.942175000000001e+28, 1.193061e+29, 1.3919045000000002e+29, 1.4316731999999999e+29, 1.4913262499999998e+29, 1.590748e+29, 1.7895915e+29, 1.9884350000000003e+29, 2.1872785e+29, 2.5849655e+29, 2.9826524999999997e+29, 3.3803395000000004e+29, 3.9768700000000005e+29, 5.9653049999999994e+29, 7.953740000000001e+29, 9.942175e+29, 1.1930609999999999e+30, 1.3919044999999998e+30, 1.5907480000000002e+30, 1.7895915e+30, 1.988435e+30, 2.1872785000000001e+30, 2.3861219999999998e+30, 2.5849655e+30, 2.7838089999999995e+30]),
                        list([7.95374e+28, 9.942175000000001e+28, 1.193061e+29, 1.3919045000000002e+29, 1.4316731999999999e+29, 1.4913262499999998e+29, 1.590748e+29, 1.7895915e+29, 1.9884350000000003e+29, 2.1872785e+29, 2.5849655e+29, 2.9826524999999997e+29, 3.3803395000000004e+29, 3.9768700000000005e+29, 5.9653049999999994e+29, 7.953740000000001e+29, 9.942175e+29, 1.1930609999999999e+30, 1.3919044999999998e+30, 1.5907480000000002e+30, 1.7895915e+30, 1.988435e+30, 2.1872785000000001e+30, 2.3861219999999998e+30, 2.5849655e+30, 2.7838089999999995e+30]),
                        list([9.942175000000001e+28, 1.193061e+29, 1.3919045000000002e+29, 1.4316731999999999e+29, 1.4913262499999998e+29, 1.590748e+29, 1.7895915e+29, 1.9884350000000003e+29, 2.1872785e+29, 2.5849655e+29, 2.9826524999999997e+29, 3.3803395000000004e+29, 3.9768700000000005e+29, 5.9653049999999994e+29, 7.953740000000001e+29, 9.942175e+29, 1.1930609999999999e+30, 1.3919044999999998e+30, 1.5907480000000002e+30, 1.7895915e+30, 1.988435e+30, 2.1872785000000001e+30, 2.3861219999999998e+30, 2.5849655e+30, 2.7838089999999995e+30]),
                        list([9.942175000000001e+28, 1.193061e+29, 1.3919045000000002e+29, 1.4316731999999999e+29, 1.4913262499999998e+29, 1.590748e+29, 1.7895915e+29, 1.9884350000000003e+29, 2.1872785e+29, 2.5849655e+29, 2.9826524999999997e+29, 3.3803395000000004e+29, 3.9768700000000005e+29, 5.9653049999999994e+29, 7.953740000000001e+29, 9.942175e+29, 1.1930609999999999e+30, 1.3919044999999998e+30, 1.5907480000000002e+30, 1.7895915e+30, 1.988435e+30, 2.1872785000000001e+30, 2.3861219999999998e+30, 2.5849655e+30, 2.7838089999999995e+30]),
                        list([9.942175000000001e+28, 1.193061e+29, 1.3919045000000002e+29, 1.4316731999999999e+29, 1.4913262499999998e+29, 1.590748e+29, 1.7895915e+29, 1.9884350000000003e+29, 2.1872785e+29, 2.5849655e+29, 2.9826524999999997e+29, 3.3803395000000004e+29, 3.9768700000000005e+29, 5.9653049999999994e+29, 7.953740000000001e+29, 9.942175e+29, 1.1930609999999999e+30, 1.3919044999999998e+30, 1.5907480000000002e+30, 1.7895915e+30, 1.988435e+30, 2.1872785000000001e+30, 2.3861219999999998e+30, 2.5849655e+30, 2.7838089999999995e+30]),
                        list([1.193061e+29, 1.3919045000000002e+29, 1.4316731999999999e+29, 1.4913262499999998e+29, 1.590748e+29, 1.7895915e+29, 1.9884350000000003e+29, 2.1872785e+29, 2.5849655e+29, 2.9826524999999997e+29, 3.3803395000000004e+29, 3.9768700000000005e+29, 5.9653049999999994e+29, 7.953740000000001e+29, 9.942175e+29, 1.1930609999999999e+30, 1.3919044999999998e+30, 1.5907480000000002e+30, 1.7895915e+30, 1.988435e+30, 2.1872785000000001e+30, 2.3861219999999998e+30, 2.5849655e+30, 2.7838089999999995e+30]),
                        list([1.193061e+29, 1.3919045000000002e+29, 1.4316731999999999e+29, 1.4913262499999998e+29, 1.590748e+29, 1.7895915e+29, 1.9884350000000003e+29, 2.1872785e+29, 2.5849655e+29, 2.9826524999999997e+29, 3.3803395000000004e+29, 3.9768700000000005e+29, 5.9653049999999994e+29, 7.953740000000001e+29, 9.942175e+29, 1.1930609999999999e+30, 1.3919044999999998e+30, 1.5907480000000002e+30, 1.7895915e+30, 1.988435e+30, 2.1872785000000001e+30, 2.3861219999999998e+30, 2.5849655e+30, 2.7838089999999995e+30]),
                        list([1.3919045000000002e+29, 1.4316731999999999e+29, 1.4913262499999998e+29, 1.590748e+29, 1.7895915e+29, 1.9884350000000003e+29, 2.1872785e+29, 2.5849655e+29, 2.9826524999999997e+29, 3.3803395000000004e+29, 3.9768700000000005e+29, 5.9653049999999994e+29, 7.953740000000001e+29, 9.942175e+29, 1.1930609999999999e+30, 1.3919044999999998e+30, 1.5907480000000002e+30, 1.7895915e+30, 1.988435e+30, 2.1872785000000001e+30, 2.3861219999999998e+30, 2.5849655e+30, 2.7838089999999995e+30]),
                        list([1.3919045000000002e+29, 1.4316731999999999e+29, 1.4913262499999998e+29, 1.590748e+29, 1.7895915e+29, 1.9884350000000003e+29, 2.1872785e+29, 2.5849655e+29, 2.9826524999999997e+29, 3.3803395000000004e+29, 3.9768700000000005e+29, 5.9653049999999994e+29, 7.953740000000001e+29, 9.942175e+29, 1.1930609999999999e+30, 1.3919044999999998e+30, 1.5907480000000002e+30, 1.7895915e+30, 1.988435e+30, 2.1872785000000001e+30, 2.3861219999999998e+30, 2.5849655e+30]),
                        list([1.3919045000000002e+29, 1.4316731999999999e+29, 1.4913262499999998e+29, 1.590748e+29, 1.7895915e+29, 1.9884350000000003e+29, 2.1872785e+29, 2.5849655e+29, 2.9826524999999997e+29, 3.3803395000000004e+29, 3.9768700000000005e+29, 5.9653049999999994e+29, 7.953740000000001e+29, 9.942175e+29, 1.1930609999999999e+30, 1.3919044999999998e+30, 1.5907480000000002e+30, 1.7895915e+30, 1.988435e+30, 2.1872785000000001e+30, 2.3861219999999998e+30]),
                        list([1.3919045000000002e+29, 1.4316731999999999e+29, 1.4913262499999998e+29, 1.590748e+29, 1.7895915e+29, 1.9884350000000003e+29, 2.1872785e+29, 2.5849655e+29, 2.9826524999999997e+29, 3.3803395000000004e+29, 3.9768700000000005e+29, 5.9653049999999994e+29, 7.953740000000001e+29, 9.942175e+29, 1.1930609999999999e+30, 1.3919044999999998e+30, 1.5907480000000002e+30, 1.7895915e+30, 1.988435e+30, 2.1872785000000001e+30]),
                        list([1.3919045000000002e+29, 1.4316731999999999e+29, 1.4913262499999998e+29, 1.590748e+29, 1.7895915e+29, 1.9884350000000003e+29, 2.1872785e+29, 2.5849655e+29, 2.9826524999999997e+29, 3.3803395000000004e+29, 3.9768700000000005e+29, 5.9653049999999994e+29, 7.953740000000001e+29, 9.942175e+29, 1.1930609999999999e+30, 1.3919044999999998e+30, 1.5907480000000002e+30, 1.7895915e+30, 1.988435e+30]),
                        list([1.3919045000000002e+29, 1.4316731999999999e+29, 1.4913262499999998e+29, 1.590748e+29, 1.7895915e+29, 1.9884350000000003e+29, 2.1872785e+29, 2.5849655e+29, 2.9826524999999997e+29, 3.3803395000000004e+29, 3.9768700000000005e+29, 5.9653049999999994e+29, 7.953740000000001e+29, 9.942175e+29, 1.1930609999999999e+30, 1.3919044999999998e+30, 1.5907480000000002e+30, 1.7895915e+30])], dtype='object')
    return (ages, t_effs, radii, masses)

def simulate_l_t(SpT):
    '''Given a spectral type (and an assumption of a uniform random age between 0.005 and 10 Gyr), 
        returns an L/T radius (in m) and mass (in kg) based on Baraffe models. 
    '''

    #Get the temperature based on the SpT
    t_eff = faherty_2016_spt_teff_relation(SpT)
    #Draw a random age from 0.005 to 10 Gyr
    age = np.random.uniform(0.005, 10)
    
    b_ages, b_teffs, b_radii, b_masses = baraffe_models()
    
    age_ind = np.where(abs(age - b_ages) == np.min(abs(age - b_ages)))[0][0]
    t_eff_arr = np.array(b_teffs[age_ind])

    t_eff_ind = np.where(abs(t_eff - t_eff_arr) == np.min(abs(t_eff - t_eff_arr)))[0][0]
    
    radius = b_radii[age_ind][t_eff_ind]
    mass = b_masses[age_ind][t_eff_ind]
    
    return mass, radius

def calc_a(period,mStar):
    G = const.G.value
    #Calculates a, the planet's semi-major axis, from its orbital period P in days and the host mass in kg. 
    a = (G*mStar*(period*86400)**2/(4*np.pi**2))**(1/3)
    return a
