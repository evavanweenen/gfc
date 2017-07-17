import gfc

sqwrange = gfc.np.arange(0.5, 3.75, 0.25)
wrange = sqwrange**2.
Krange = range(4, 21)[13:]
for w in wrange:
    f_w = "tgas/2MASS/XD_K_w/w_{0}".format(w)
    if not gfc.isdir(f_w):
        gfc.mkdir(f_w)
    for K in Krange:
        f = "tgas/2MASS/XD_K_w/w_{0}/K_{1}".format(w, K)
        if not gfc.isdir(f):
            gfc.mkdir(f)

t = gfc.io.load_table_with_separate_arrays(saveto_folder="tgas/2MASS/results/")

warr = gfc.XD_arr(t, "w1", "w2", "w3")
wcov = gfc.XD_arr(t, "S") ; wcov[:, 0, 0] = 1e15
proj = gfc.XD_arr(t, "R")

initial_amps = [23.92, 22.78, 22.62, 8.91, 8.51, 6.74, 3.98, 1.74, 0.72, 0.09, 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.] # %

initial_means = [(5.54, -6.97, -9.26),
                 (-22.72, -10.24, -7.3),
                 (-12.92, -29.57, -7.6),
                 (-19.43, -32.9, -4.93),
                 (9.2, 3.89, -5.99),
                 (-17.79, -22.69, -4.51),
                 (-9.07, -20.50, -4.88),
                 (-40.07, -18.92, 0.64),
                 (-28.28, -105.62, 2.87),
                 (2.08, -103.07, -8.2),
                 (1., 2., 3.),
                 (10., 5., -10.),
                 (20., 10., -35.),
                 (30., 15., -50.),
                 (40., 20., 20.),
                 (50., -30., 30.),
                 (60., -20., -20.),
                 (-10., -10., 15.),
                 (-20., -5., 60.),
                 (-30., 45., 0.)]

initial_covs = [(700., -111., -60., 200., 25., 145.),
                (243, 70, 10, 48, 10, 40),
                (1836, 55, -60, 670, -30, 540),
                (350, 165, 110, 230, 110, 135),
                (80, -30, -20, 25, 9, 50),
                (70, -30, -45, 18, 17, 60),
                (9, -7, -6, 13, 8, 14),
                (30, 0.2, 13, 0.6, 0.2, 6),
                (4600, -2400, -330, 4500, -125, 3500),
                (1, 0.5, 0.3, 4, 0.5, 3),
                (19., 0., 0., 10., 0., 10.),
                (18., 0., 0., 11., 0., 10.),
                (17., 0., 0., 12., 0., 10.),
                (16., 0., 0., 13., 0., 10.),
                (15., 0., 0., 14., 0., 10.),
                (14., 0., 0., 15., 0., 10.),
                (13., 0., 0., 16., 0., 10.),
                (12., 0., 0., 17., 0., 10.),
                (11., 0., 0., 18., 0., 10.),
                (10., 0., 0., 19., 0., 10.),]
initial_covs = [gfc.covariant_array(cov) for cov in initial_covs]

for K in Krange:
    for w in wrange:
        f = "tgas/2MASS/XD_K_w/w_{0}/K_{1}".format(w, K)
        print "K = {0} ; w = {1}".format(K, w),  ; g.io.flush()
        amps_xd, means_xd, covs_xd, L = g.XD(warr, wcov, initial_amps[:K], initial_means[:K], initial_covs[:K], projection = proj, w = w)
        print >> open(f+"/L", 'w'), L
        gfc.io.save_PDFs(amps_xd, means_xd, covs_xd, f)
        print "done!"
