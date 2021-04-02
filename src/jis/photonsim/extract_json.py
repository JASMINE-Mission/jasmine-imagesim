import numpy as np
import os
import json
import dataclasses
from jis.pixsim import makeflat as mf
from jis.pixsim import readflat as rf
from jis.photonsim import aperture

def extsp(sp):
    """
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
    This function extracts detector properties from an input json file
    and returns a detector class object which has the properties..

    Args:
        det_json_filename (str): Filename of the detector json file.
        spixdim    (array-like): Dimensions of the intrapixel pattern.
                                 (Default: [32, 32])

    Returns:
        detector (detector): Detector object which has the extracted properties.

    """

    @dataclasses.dataclass(frozen=True)
    class Persistence:
        """
        This class contains the detector persistence parameters.

        Attributes:
            tau (ndarray): Detrapping timescales in seconds.
            rho (ndarray): Trapping fractions.
        """
        tau : np.ndarray
        rho : np.ndarray


    @dataclasses.dataclass(frozen=True)
    class ReadParams:
        """
        This class contains the detector readout parameters.

        Attributes:
            fsmpl      (float): Sampling frequency in Hz.
            tsmpl (float)          : Sampling time in sec (1/fsmpl).
            ncol_ch (int)          : Number of columns in one ch.
            nrow_ch (int)          : Number of rows in one ch.
            npix_pre (int)         : Number of pixels before reading each row.
            npix_post (int)        : Number of pixels after reading each row.
            t_overhead (float)     : Overhead time between reset and the 1st read in sec.
            npix_read_per_row (int): Number of pixels in a row including pre/post pixels.
        """
        fsmpl            : float = None
        tsmpl            : float = dataclasses.field(init=False)
        ncol_ch          : int   = None
        nrow_ch          : int   = None
        npix_pre         : int   = None
        npix_post        : int   = None
        npix_read_per_row: int = dataclasses.field(init=False)
        t_overhead       : float = None
        t_scan           : float = dataclasses.field(init=False)

        def __post_init__(self):
            tsmpl = 1.0/self.fsmpl if self.fsmpl else None
            npix_read_per_row = self.ncol_ch + self.npix_pre + self.npix_post
            t_scan = self.t_overhead + tsmpl * npix_read_per_row * self.nrow_ch
            object.__setattr__(self, 'tsmpl', tsmpl)
            object.__setattr__(self, 'npix_read_per_row', npix_read_per_row)
            object.__setattr__(self, 't_scan', t_scan)


    @dataclasses.dataclass(frozen=True)
    class QuantumEfficiency:
        """
        This class defines the detector quantum efficiency.

        Attributes:
            wl (ndarray) : Wavelengthes in um.
            val (ndarray): Quantum efficiencies.
        """
        wl : np.ndarray
        val: np.ndarray


    @dataclasses.dataclass(frozen=True)
    class FlatFrame:
        """
        This class contains the inter- and intra-pixel flat frames.
        """
        intrapix   : np.ndarray
        interpix   : np.ndarray


    @dataclasses.dataclass(frozen=True)
    class Detector:
        """
        This is a class to describe the detector properties.

        Attributes:
            npix (int)               : Number of pixels on a side.
            idark (float)            : Dark current including stray light (e/s/pix).
            intrapix (ndarray)       : Intrapixel pattern.
            interpix (ndarray)       : Interpixel pattern.
            readnoise (float)        : Readnoise (e/read).
            fullwell (float)         : Full-well in electrons.
            gain (float)             : Conversion gain in e/adu.
            pixsize (float)          : Pixel size in um.
            qe (QuantumEfficiency)   : Quantum efficiency.
            persistence (Persistence): Persistence parameters consists of tau and rho.
            readparams (ReadParams)  : Detector readout parameters.
        """
        npix       : int
        idark      : float
        readnoise  : float
        fullwell   : float
        gain       : float
        pixsize    : float
        flat       : FlatFrame
        qe         : QuantumEfficiency
        persistence: Persistence
        readparams : ReadParams


    with open(det_json_filename, "r") as fp:
        js = json.load(fp)

    npix = js['Npix']['val']
    tau, rho = extPersistenceParams(js)
    wl_qe, val_qe = extQE(js)
    interpix_sigma, intrax, intray = extFlatInfo(js)

    detector = Detector(
        npix       = npix,
        idark      = extIdark(js),
        readnoise  = js['readnoise']['val'],
        fullwell   = js['Fullwell']['val'],
        gain       = js['gain']['val'],
        pixsize    = js['pixsize']['val'],
        flat       = FlatFrame(
            interpix = mf.gaussian_flat(Nside=npix, sigma=interpix_sigma),
            intrapix = rf.read_intrapix(intrax, intray, spixdim)),
        qe         = QuantumEfficiency(wl = wl_qe, val = val_qe),
        persistence= Persistence(tau = tau, rho = rho),
        readparams = ReadParams(
            fsmpl      = js['readparams']['fsmpl']['val'],
            ncol_ch    = js['readparams']['ncol_ch']['val'],
            nrow_ch    = js['readparams']['nrow_ch']['val'],
            npix_pre   = js['readparams']['npix_pre']['val'],
            npix_post  = js['readparams']['npix_post']['val'],
            t_overhead = js['readparams']['t_overhead']['val']))

    return detector


def extIdark(det):
    """
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
    This function extracts persistence parameters.

    Args:
        det (dict): Decoded json data of the detector information.

    Returns:
        tau (ndarray): Detrapping time constants in sec.
        rho (ndarray): Trapping fractions.

    """

    tau = np.array(det['persistence']['tau']['val'])
    rho = np.array(det['persistence']['rho']['val'])

    return tau, rho


def extFlatInfo(det):
    """
    This function extracts information about detector flat.

    Args:
        det (dict): Decoded json data of the detector information.

    Returns:
        interpix (float): Stddev of the interpixel flat.
        intrax   (str)  : Filename of the intrapix pattern in x direction.
        intray   (str)  : Filename of the intrapix pattern in y direction.

    """

    interpix = det['interpix']['stddev']['val']
    intrax   = det['intrapix']['file_x']['val']
    intray   = det['intrapix']['file_y']['val']

    return interpix, intrax, intray


def extQE(det):
    """
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
    wl = np.array(det['QE']['wavelength']['val'])
    qe = np.array(det['QE']['qe_value']['val'])

    return wl, qe


def extsrc(src):
    """
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
    This method creates an instance of ControlParams,
    which contains control parameters written in the
    json file for the control parameters.

    Args:
        json_filename (string): Input json filename.

    Returns:
        ControlParams: Created ControlParams object.

    """

    @dataclasses.dataclass(frozen=True)
    class EffectSelector:
        """
        This class handles the switch to enable/disable the effects.

        Attributes:
            ace (bool): incorporate the Attitude Control Error.
            wfe (bool): incorpolate the Wave Front Error.
            flat_interpix (bool): incorporate the Interpix Flat Error.
            flat_intrapix (bool): incorporate the Intrapix Flat Error.
            psf (bool): simulate a physically realistic Point Spread Function.
            rolling_shutter (bool): incorpolate the time-delay due to rolling shutter.
        """
        ace: bool
        wfe: bool
        flat_interpix: bool
        flat_intrapix: bool
        psf: bool
        rolling_shutter: bool


    @dataclasses.dataclass(frozen=True)
    class ControlParams:
        """
        This is a class to handle the control parameters.

        Attributes:
            wfe_control (dict): Parameters related to the wfe calculation.
            M_parameter (int): The M parameter which determine the cell scale of the psf.
                               The fp-cell scale will be (1/M) x 10^-3 rad/fp-cell.
            ace_control (dict): Parameters related to the ace calculation.
            nplate (int): Number of plates that make up a small frame.
            effect (EffectSelector): Flags to enable/disable components.

        """
        wfe_control: dict
        M_parameter: int
        ace_control: dict
        nplate     : int
        effect     : EffectSelector


    with open(json_filename, "r") as fp:
        js = json.load(fp)

        wfe = {}
        wfe_control = js.get('WFEcontrol')
        wfe_items = [
            'zernike_nmax', 'zernike_even', 'zernike_odd', 'reference_wl']
        if wfe_control:
            for item in wfe_items:
                wfe[item] = wfe_control[item]['val']

        M = js['M']['val']

        ace = {}
        ace_control = js.get('ACEcontrol')
        ace_items = ['dtace', 'tace', 'acex_std', 'acey_std']
        if ace_control:
            for item in ace_items:
                ace[item] = ace_control[item]['val']
            ace['nace'] = int(ace['tace']/ace['dtace'])+1
            ace['tace'] = ace['dtace']*ace['nace']

        nplate = js['Nplate']['val']

        effect_obj = js.get('effect')
        effect = {}
        for field in dataclasses.fields(EffectSelector):
            try:
                item = effect_obj.get(field.name)
                effect[field.name] = item.get('val', False)
            except:
                effect[field.name] = False
        effect = EffectSelector(**effect)

    control_params = ControlParams(
        wfe_control=wfe, M_parameter=M, ace_control=ace,
        nplate=nplate, effect=effect)

    return control_params


def mkTel(json_filename):
    """
    This method creates a telescope object.

    Args:
        json_filename: Input json filename.

    Returns:
        telescope: Created telescope object.

    """

    class telescope:
        """
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
