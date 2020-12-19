import numpy as np
import os 
import json
from jis.pixsim import makeflat as mf
from jis.pixsim import readflat as rf
from jis.photonsim import aperture

def extsp(sp):
    """
    Summary:
        This function extracts spectral information 
        from a decoded spectral data, 'sp'.

    Args:
        sp (dict): Decoded json data of the spectral information.

    Returns:
        k    (int)    : Number of output data points.
        WL   (ndarray): Wavelength data (um).
        NP   (ndarray): Number of detected photon (electron) flux (e-/s/m^2/um).
        Ntot (float)  : Total detected photon (electron) rate (e-/s/m^2).

    Example:
        import json
        from jis.photonsim.extract_json import extsp

        fp = open("***.json")
        sp = json.load(fp) 
        k, WL, NP, Ntot = extsp(sp)

    """

    k  = len(sp['WLdef'])-1

    WL = np.empty(k)
    for i in range(k):
        WL[i] = sp['WLdef']['v{:02d}'.format(i)]

    NP = np.empty(k)
    for i in range(k):
        NP[i] = sp['SPR']['v{:02d}'.format(i)]

    Ntot = sp['Ntot']['val']

    return k, WL, NP, Ntot


def exttel(tel):
    """
    Summary:
        This function extracts telescope information
        from a decoded telescope data, 'tel'.

    Args:
        tel (dict): Decoded json data of the telescope information.

    Returns:
        k         (int)    : Number of output data points.
        WLdefined (ndarray): Wavelength data (um).
        EPdefined (ndarray): Optical efficiency.
        WLshort   (float)  : Shortest wavelength (um).
        WLlong    (float)  : Longest wavelength (um).        

    Example:
        import json
        from jis.photonsim.extract_json import exttel

        fp  = open("***.json")
        tel = json.load(fp) 
        k, WL, EP, WLmin, WLmax = exttel(tel)

    """

    # from Eopt in tel, prepare arrays for WL and E
    k = int((len(tel['Eopt'])-1)/2 ) # number of definition points

    WLdefined = np.empty(k) # Wavelength
    EPdefined = np.empty(k) # Efficiency of optics
    for i in range(k):
        WLdefined[i] = tel['Eopt']['W{:02d}'.format(i)]
        EPdefined[i] = tel['Eopt']['V{:02d}'.format(i)]

    WLshort = WLdefined[0]
    WLlong =  WLdefined[k-1]

    return k, WLdefined, EPdefined, WLshort, WLlong


def mkDet(det_json_filename, spixdim=[32, 32]):
    """
    Summary:
        This function extracts detector properties from an input json file
        and returns a detector class object which has the properties..

    Args:
        det_json_filename (str): Filename of the detector json file.
        spixdim    (array-like): Dimensions of the intrapixel pattern.
                                 (Default: [32, 32])

    Returns:
        detector (detector): Detector object which has the extracted properties.

    """

    class detector:
        def __init__(self, npix=None, idark=None, intrapix=None, interpix=None,\
                     tau=None, rho=None, readnoise=None, fullwell=None,\
                     gain=None, readparams=None, qe=None, pixsize=None):
            """
            Summary:
                This is a class to describe the detector properties.

            Attributes:
                npix     (int)    : Number of pixels on a side.
                idark    (float)  : Dark current including stray light (e/s/pix).
                intrapix (ndarray): Intrapixel pattern.
                interpix (ndarray): Interpixel pattern.
                persistence (dict): Persistence parameters consists of tau and rho.
                                    tau is the detrapping timescales in sec.
                                    rho is the trapping fractions.
                readnoise  (float): Readnoise (e/read).
                fullwell   (float): Full-well in electrons.
                gain       (float): Conversion gain in e/adu.
                fsmpl      (float): Sampling frequency in Hz.
                tsmpl      (float): Sampling time in sec (1/fsmpl).
                ncol_ch    (int)  : Num. of col. in one ch.
                nrow_ch    (int)  : Num. of row in one ch.
                npix_pre   (int)  : Npix before reading each row.
                npix_post  (int)  : Npix after reading each row.
                t_overhead (float): Overhead time between reset and the 1st read in sec.
                qe         (dict) : Quantum efficiency (wl: wavelength in um; val: qe values).
                pixelsize  (float): Pixel size in um.

            """

            self.npix  = npix
            self.idark = idark
            self.intrapix = intrapix
            self.interpix = interpix
            self.persistence = {'tau': np.array(tau), 'rho': np.array(rho)}
            self.readnoise = readnoise
            self.fullwell = fullwell
            self.gain = gain
            self.qe = qe
            self.pixsize=pixsize
            if readparams is not None:
                self.fsmpl      = readparams['fsmpl']['val']
                self.tsmpl      = 1./self.fsmpl
                self.ncol_ch    = readparams['ncol_ch']['val']
                self.nrow_ch    = readparams['nrow_ch']['val']
                self.npix_pre   = readparams['npix_pre']['val']
                self.npix_post  = readparams['npix_post']['val']
                self.t_overhead = readparams['t_overhead']['val']
            else:
                self.fsmpl      = None
                self.tsmpl      = None
                self.ncol_ch    = None
                self.nrow_ch    = None
                self.npix_pre   = None
                self.npix_pos   = None
                self.t_overhead = None


    with open(det_json_filename, "r") as fp:
        js = json.load(fp)

        idark    = extIdark(js)
        tau, rho = extPersistenceParams(js)
        interpix_sigma, intradir, intrax, intray = extFlatInfo(js)

        npix       = js['Npix']['val']
        readnoise  = js['readnoise']['val']
        readparams = js['readparams']
        fullwell   = js['Fullwell']['val']
        gain       = js['gain']['val']
        pixsize    = js['pixsize']['val']

        wl_qe, val_qe = extQE(js)
    fp.close()

    qe = {'wl': np.array(wl_qe), 'val': np.array(val_qe)}

    interpix = mf.gaussian_flat(Nside=npix, sigma=interpix_sigma)
    intrapix = rf.read_intrapix(intrax, intray, spixdim, intradir)

    detector = detector(npix=npix, idark=idark, \
                        interpix=interpix, intrapix=intrapix, \
                        tau=tau, rho=rho, readnoise=readnoise,\
                        fullwell=fullwell, gain=gain,\
                        readparams=readparams, qe=qe, pixsize=pixsize)

    return detector


def extIdark(det):
    """
    Summary:
        This function extracts dark current including stray light.

    Args:
        det (dict): Decoded json data of the detector information.

    Returns:
        idark (float): Dark current including stray light in e/sec/pix.

    """

    idark = det['Idark']['val']

    return idark


def extPersistenceParams(det):
    """
    Summary:
        This function extracts persistence parameters.

    Args:
        det (dict): Decoded json data of the detector information.

    Returns:
        tau (ndarray): Detrapping time constants in sec.
        rho (ndarray): Trapping fractions.

    """

    tau = det['persistence']['tau']
    rho = det['persistence']['rho']

    return tau, rho


def extFlatInfo(det):
    """
    Summary:
        This function extracts information about detector flat.

    Args:
        det (dict): Decoded json data of the detector information.

    Returns:
        interpix (float): Stddev of the interpixel flat.
        intradir (str)  : Name of the directory containing intrapix pattern data.
        intrax   (str)  : Filename of the intrapix pattern in x direction.
        intray   (str)  : Filename of the intrapix pattern in y direction.

    """

    interpix = det['interpix']['val']

    intradir = det['intrapix']['dirname']
    intrax   = det['intrapix']['filex']
    intray   = det['intrapix']['filey']

    return interpix, intradir, intrax, intray


def extQE(det):
    """
    Summary:
        This function extracts quantum efficiency information
        from a decoded detector data, 'det'.

    Args:
        det (dict): Decoded json data of the detector information.

    Returns:
        WLdet (ndarray): Wavelength data (um).
        QEdet (ndarray): Quantum efficiency of the detector.

    Example:
        import json
        from jis.photonsim.extract_json import extQE

        fp  = open("***.json")
        det = json.load(fp) 
        WL, QE = extQE(det)

    """ 

    # read QE from det
    k = int((len(det['QE'])-1)/2 ) # number of definition points

    WLdet = np.empty(k) # Wave length
    QEdet = np.empty(k) # Efficiency of optics
    for i in range(k):
        WLdet[i] = det['QE']['W{:02d}'.format(i)]
        QEdet[i] = det['QE']['V{:02d}'.format(i)]

    return WLdet, QEdet


def extsrc(src):
    """
    Summary:
        This function extracts object source information
        from a decoded object source data, 'src'.

    Args:
        src (dict): Decoded json data of the object source information.

    Returns:
        Rv  (float): Extinction factor Rv(=Av/E(B-V)).
        JH  (float): Color excess (E(J-H)=AJ-AH).
        alp (float): Interpolation factor to define the Hw-band mag
                         Hw-mag = alp * J-mag + (1-alp) * H-mag

    Example:
        import json
        from jis.photonsim.extract_json import extsrc

        fp  = open("***.json")
        src = json(fp) 
        Rv, JH, alp = extsrc(src) 

    """
    
    Rv  = src['Rv']['val']
    JH  = src['J-H']['val']
    alp = src['Hwband']['val']

    return Rv, JH, alp


def mkControlParams(json_filename):
    """
    Summary:
        This method creates a control_params object,
        which contains control parameters written in the 
        json file for the control parameters.

    Args:
        json_filename (string): Input json filename.

    Returns:
        control_params (control_params): Created control_params object.

    """


    class control_params:
        """
        Summary:
            This is a class to handle the control parameters.

        Attributes:
            wfe_control (dict): Parameters related to the wfe calculation.
            M_parameter (int) : The M parameter which determine the cell scale of the psf.
                                The fp-cell scale will be (1/M) x 10^-3 rad/fp-cell.
            ace_control (dict): Parameters related to the ace calculation.

        """
          
        def __init__(self, wfe=None, M=None, ace=None):
            self.wfe_control = wfe
            self.M_parameter = M
            self.ace_control = ace


    with open(json_filename, "r") as fp:
        js = json.load(fp)

        wfe = {}
        for item in ['zernike_nmax', 'zernike_even', 'zernike_odd', 'reference_wl']:
            wfe[item] = js['WFEcontrol'][item]['val']

        M = js['M']['val']

        ace = {}
        for item in ['dtace', 'tace']:
            ace[item] = js['ACEcontrol'][item]['val'] 
        ace['nace'] = int(ace['tace']/ace['dtace'])+1
        ace['tace'] = ace['dtace']*ace['nace']
    fp.close()

    control_params = control_params(wfe=wfe, M=M, ace=ace)

    return control_params


def mkTel(json_filename):
    """
    Summary:
        This method creates a telescope object.

    Args:
        json_filename: Input json filename.

    Returns:
        telescope: Created telescope object.

    """

    class telescope:
        """
        Summary:
            This is a class to handle telescope parameters.

        Attributs:
            epd              (float)  : Exit pupil diameter (mm).
            aperture         (ndarray): Aperture pattern array (ap-cell size = 1 mm).
            cobs             (float)  : Obscuration ratio (R(obscuration)*2/EPD).
            efl              (float)  : Effective focal length (mm).
            spider_type      (string) : Spider type.
            spider_thickness (float)  : Spider thickness (mm).
            total_area       (float)  : Telescope total area (m^2).
            opt_efficiency   (dict)   : Optical efficiency.

        """
           
        def __init__(self, epd=None, aperture=None,\
                     cobs=None, spider=None, total_area=None,\
                     opt_efficiency=None, efl=None):
            self.epd      = epd                         # Exit pupil diameter (mm).
            self.aperture = aperture                    # Aperture pattern.
            self.cobs     = cobs                        # Obscuration ratio (Robs*2/EPD).
            self.efl      = efl                         # Effective focal length (mm).
            self.spider_type      = spider['type']      # Spider tape.
            self.spider_thickness = spider['thickness'] # Spider thickness in mm.
            self.total_area       = total_area          # Total area in m^2.
            self.opt_efficiency   = opt_efficiency      # Optics efficiency.


    with open(json_filename, "r") as fp:
        js = json.load(fp)

        epd  = js['EPD']['val']     # Exit pupil diameter in mm.
        cobs = js['Cobs']['val']    # Obscuration ratio.
        efl  = js['EFL']['val']     # Effective focal length in mm.
        r_obscuration = epd/2.*cobs # Obscuration radius in mm.
        spider_params = {'type':js['Stype']['val'], 'thickness':js['Stype']['thick']}
        n_opteff, wl_opteff, opteff, wl_opteff_short, wl_opteff_long = exttel(js)
    fp.close()

    opt_efficiency = {'wl': np.array(wl_opteff), 'val': np.array(opteff)}

    ap_data    = None
    total_area = None
    if spider_params['type'] == 'tripod':
        n_apcell = int(epd+4)   # Assuming ap-cell scale to be 1mm/ap-cell.
                                # Set the aperture pattern size to be 2mm larger than D.

        if n_apcell%2 == 1: # n_apcell should be even.
            n_apcell = n_apcell + 1

        ap_data, total_area_mm2 = aperture.calc_aperture(n_apcell, epd, r_obscuration,\
                                                         spider_params['thickness'])

    telescope = telescope(epd=epd, aperture=ap_data, cobs=cobs,\
                          spider=spider_params, total_area=total_area_mm2*1.e-6,\
                          opt_efficiency=opt_efficiency, efl=efl)
    # total_area is in m^2.

    return telescope

def mkVar(json_filename):
    """
    Summary: This is dummy to initialize variability class, to match mkDet etc.
    """            
    _var = variability(json_filename)                    
    return _var

class variability():
    """
    Summary: stellar variability class.
    """
    def __init__(self,json_filename):
        self.var = True
        self.read_json(json_filename)

    def read_json(self,json_filename):
        """
        Summary: this function is json i/o for stellar variability
        """
        with open(json_filename, "r") as f:
            var_params = json.load(f)            
            f.close()
            self.Nvar=len(var_params)
            self.plate=[]
            self.star=[]
            self.vartype=[]
            self.varfile=[]
            self.dirname=[]
            for i in range(0,self.Nvar):
                var=var_params[str(i)]
                self.plate.append(var["plate"])
                self.star.append(var["star"])
                self.vartype.append(var["vartype"])
                self.varfile.append(var["varfile"])
                self.dirname.append(var["dirname"])
    
    def read_var(self,t_day,star_index):
        """
        Summary
        -------        
        read variability for a given star index

        Parameters
        ----------
        t_day : time array in the unit of day
        star_index : star index

        Returns
        -------
        sw: if True there is variability, else no variability.
        injlc: variability
        b : impact parameter for vartype = planet
       
        """
        from jis.pixsim import transitmodel 
        for i in range(0,self.Nvar):
            if  int(star_index) == int(self.star[i]):
                if self.vartype[i] == "planet":
                    injlc, b=transitmodel.gentransit_json(t_day,os.path.join(self.dirname[i],self.varfile[i]))
                    sw=True
                    return sw,injlc,b
                else:
                    sys.exit("No valid vartype")

        sw=False
        return sw,None,None

            
