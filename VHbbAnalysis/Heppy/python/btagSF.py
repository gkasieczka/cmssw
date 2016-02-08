import ROOT
import os

csvpath = os.environ['CMSSW_BASE']+"/src/VHbbAnalysis/Heppy/data/csv"
ROOT.gSystem.Load(csvpath+'/BTagCalibrationStandalone.so')
calib = ROOT.BTagCalibration("csvv2", csvpath+"/CSVv2.csv")

for wp in [ [0, "L"],[1, "M"], [2,"T"] ]:
    print "CSV"+wp[1]+":"
    for syst in ["central", "up", "down"]:
        print "\t"+syst+":"
        csv_calib_bc = ROOT.BTagCalibrationReader(calib, wp[0], "mujets", syst)
        csv_calib_l = ROOT.BTagCalibrationReader(calib, wp[0], "comb", syst)
        weight = lambda pt, eta, fl, csv_calib_bc=csv_calib_bc, csv_calib_l=csv_calib_l : (csv_calib_bc.eval( -fl+5 ,eta, pt) if (abs(eta)<2.4 and pt<670. and pt>30.) else 1.0) if fl>=4 else (csv_calib_l.eval(2,eta, pt) if (abs(eta)<2.4 and pt<1000. and pt>20.) else 1.0)
        print "\t\tB:", weight(50., 0.5, 5)
        print "\t\tC:", weight(50., 0.5, 4)
        print "\t\tL:", weight(50., 0.5, 0)



'''
reader_L = ROOT.BTagCalibrationReader(calib, 0, "mujets", "central")  # 0 is for loose op
reader_M = ROOT.BTagCalibrationReader(calib, 1, "mujets", "central")  # 0 is for loose op
reader_T = ROOT.BTagCalibrationReader(calib, 2, "mujets", "central")  # 0 is for loose op

reader2_L = ROOT.BTagCalibrationReader(calib, 0, "comb", "central")  # 0 is for loose op
reader2_M = ROOT.BTagCalibrationReader(calib, 1, "comb", "central")  # 0 is for loose op
reader2_T = ROOT.BTagCalibrationReader(calib, 2, "comb", "central")  # 0 is for loose op


# in your event loop
for fl in [0,1,2]:
    if(fl==2):
        print str(fl)+", loose: " , reader2_L.eval(fl, 1.2, 300000.)  # jet flavor, eta, pt
        print str(fl)+", medium: " , reader2_M.eval(fl, 999, 30.)  # jet flavor, eta, pt
        print str(fl)+", tight: " , reader2_T.eval(fl, 1.2, 30000.)  # jet flavor, eta, pt
    else:
        print str(fl)+", loose: " , reader_L.eval(fl, 999, 70.)  # jet flavor, eta, pt
        print str(fl)+", medium: " , reader_M.eval(fl, 1.2, 90909090.)  # jet flavor, eta, pt
        print str(fl)+", tight: " , reader_T.eval(fl, 1.2, 70.)  # jet flavor, eta, pt
'''
