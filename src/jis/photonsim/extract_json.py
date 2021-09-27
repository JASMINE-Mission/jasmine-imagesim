import numpy as np
import os
import json
import dataclasses
from typing import List
from jis.pixsim import makeflat as mf
from jis.pixsim import readflat as rf
from jis.photonsim import aperture

def extSp(sp):
    """
    This function extracts spectral information
    from a decoded spectral data, 'sp'.

    Args:
        sp (dict): Decoded json data of the spectral information.

    Returns:
        wavelength (ndarray): Wavelength data (um).
        spectrum   (ndarray): Number of detected photon (electron) flux (e-/s/m^2/um).
        Ntotal     (float)  : Total detected photon (electron) rate (e-/s/m^2).

    Example:
        import json
        from jis.photonsim.extract_json import extSp

        fp = open("***.json")
        sp = json.load(fp)
        WL, NP, Ntot = extSp(sp)

    """

    wavelength = np.array(sp['wavelength']['val'])
    spectrum   = np.array(sp['spectrum']['val'])
    Ntotal = sp['Ntot']['val']

    return wavelength, spectrum, Ntotal


def extTel(tel):
    """
    This function extracts telescope optical efficiency
    from a decoded telescope data, 'tel'.

    Args:
        tel (dict): Decoded json data of the telescope information.

    Returns:
        wavelength (ndarray): Wavelength data (um).
        efficiency (ndarray): Optical efficiency.

    Example:
        import json
        from jis.photonsim.extract_json import extTel

        fp  = open("***.json")
        tel = json.load(fp)
        WL, EP = extTel(tel)

    """

    wavelength = np.array(tel['Eopt']['wavelength']['val'])
    efficiency = np.array(tel['Eopt']['efficiency']['val'])

    return wavelength, efficiency


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

    Attributes:
        intrapix (ndarray)       : Intrapixel pattern.
        interpix (ndarray)       : Interpixel pattern.
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
        readnoise (float)        : Readnoise (e/read).
        fullwell (float)         : Full-well in electrons.
        gain (float)             : Conversion gain in e/adu.
        pixsize (float)          : Pixel size in um.
        nmargin (int)            : Number of margin pixels for simpix calculation.
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
    nmargin    : np.ndarray
    flat       : FlatFrame
    qe         : QuantumEfficiency
    persistence: Persistence
    readparams : ReadParams


    @classmethod
    def from_json(cls, det_json_filename):
        """
        This function extracts detector properties from an input json file
        and returns a detector class object which has the properties..

        Args:
            det_json_filename (str): Filename of the detector json file.

        Returns:
            detector (detector): Detector object which has the extracted properties.
        """

        with open(det_json_filename, "r") as fp:
            js = json.load(fp)

        npix = js['Npix']['val']
        spixdim = np.array(js['spixdim']['val'])
        tau = np.array(js['persistence']['tau']['val'])
        rho = np.array(js['persistence']['rho']['val'])
        wl_qe, val_qe = extQE(js)
        interpix_sigma = js['interpix']['stddev']['val']
        intrapix_filex = js['intrapix']['file_x']['val']
        intrapix_filey = js['intrapix']['file_y']['val']

        detector = Detector(
            npix       = npix,
            idark      = js['Idark']['val'],
            readnoise  = js['readnoise']['val'],
            fullwell   = js['Fullwell']['val'],
            gain       = js['gain']['val'],
            pixsize    = js['pixsize']['val'],
            nmargin    = js['Nmargin']['val'],
            flat       = FlatFrame(
                interpix = mf.gaussian_flat(Nside=npix, sigma=interpix_sigma),
                intrapix = rf.read_intrapix(intrapix_filex, intrapix_filey, spixdim)),
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


def extSrc(src):
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
        from jis.photonsim.extract_json import extSrc

        fp  = open("***.json")
        src = json(fp)
        Rv, JH, alp = extSrc(src)

    """

    Rv  = src['Rv']['val']
    JH  = src['J-H']['val']
    alp = src['Hwband']['val']

    return Rv, JH, alp


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
    """
    ace: bool
    wfe: bool
    flat_interpix: bool
    flat_intrapix: bool
    psf: bool


@dataclasses.dataclass(frozen=True)
class ControlParams:
    """
    This is a class to handle the control parameters.

    Attributes:
        wfe_control (dict): Parameters related to the wfe calculation.
        M_parameter (int): The M parameter which determine the cell scale of the psf.
                           The fp-cell scale will be (1/M) x 10^-3 rad/fp-cell.
        fN_parameter (int): Number of fp-cells of the psf data calculated by calc_psf.
        ace_control (dict): Parameters related to the ace calculation.
        nplate (int): Number of plates that make up a small frame.
        tplate (float): Exposure time of each plate in second.
        Rv (float): Total-to-selective extinction of the field.
        JH (float): Color excess E(J-h) for each source.
        alpha (float): Hw-band interpolation factor.
        effect (EffectSelector): Flags to enable/disable components.
    """
    wfe_control: dict
    M_parameter: int
    fN_parameter: int
    ace_control: dict
    nplate     : int
    tplate     : float
    Rv         : float
    JH         : float
    alpha      : float
    effect     : EffectSelector

    @classmethod
    def from_json(cls, ctl_json_filename):
        """
        This method creates an instance of ControlParams,
        which contains control parameters written in the
        json file for the control parameters.

        Args:
            ctl_json_filename (string): Input json filename.

        Returns:
            ControlParams: Created ControlParams object.
        """

        with open(ctl_json_filename, "r") as fp:
            js = json.load(fp)

        wfe = {}
        wfe_control = js.get('WFEcontrol')
        wfe_items = [
            'zernike_nmax', 'zernike_even', 'zernike_odd', 'reference_wl']
        if wfe_control:
            for item in wfe_items:
                wfe[item] = wfe_control[item]['val']

        M = js['M']['val']
        fN = js['fN']['val']

        ace = {}
        ace_control = js.get('ACEcontrol')
        ace_items = ['dtace', 'tace', 'acex_std', 'acey_std']
        if ace_control:
            for item in ace_items:
                ace[item] = ace_control[item]['val']
            ace['nace'] = int(ace['tace']/ace['dtace'])+1
            ace['tace'] = ace['dtace']*ace['nace']

        nplate = js['Nplate']['val']
        tplate = js['tplate']['val']
        Rv = js['Rv']['val']
        JH = js['J-H']['val']
        alpha = js['alpha']['val']

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
            wfe_control=wfe, M_parameter=M, fN_parameter=fN, ace_control=ace, nplate=nplate,
            tplate=tplate, Rv=Rv, JH=JH, alpha=alpha, effect=effect)

        return control_params



@dataclasses.dataclass(frozen=True)
class SpiderInfo:
    """
    This class defines the spider type and thickness.

    Attributes:
        type      (str)  : Spider type.
        thickness (float): Spider thickness (mm).

    """
    type: str
    thickness: float


@dataclasses.dataclass(frozen=True)
class OpticalEfficiency:
    """
    This class defines the optical efficiency of the telescope.

    Attributes:
        wavelength (ndarray): Wavelengths in micron.
        efficiency (ndarray): Efficiencies.
    """
    wavelength: np.ndarray
    efficiency: np.ndarray


@dataclasses.dataclass(frozen=True)
class Telescope:
    """
    This is a class to handle telescope parameters.

    Attributs:
        epd            (float)  : Exit pupil diameter (mm).
        aperture       (ndarray): Aperture pattern array (ap-cell size = 1 mm).
        cobs           (float)  : Obscuration ratio (R(obscuration)*2/EPD).
        efl            (float)  : Effective focal length (mm).
        total_area     (float)  : Telescope total area (m^2).
        spider         (Spider) : Spider type.
        opt_efficiency (OpticalEfficiency) : Optical efficiency.
    """

    epd           : float
    aperture      : np.ndarray
    cobs          : float
    efl           : float
    total_area    : float
    spider        : SpiderInfo
    opt_efficiency: OpticalEfficiency


    @classmethod
    def from_json(cls, tel_json_filename):
        """
        This method creates a telescope object.

        Args:
            tel_json_filename: Input json filename.

        Returns:
            Telescope: Created telescope object.

        """

        with open(tel_json_filename, "r") as fp:
            js = json.load(fp)

        epd  = js['EPD']['val']     # Exit pupil diameter in mm.
        cobs = js['Cobs']['val']    # Obscuration ratio.
        efl  = js['EFL']['val']     # Effective focal length in mm.
        r_obscuration = epd/2.*cobs # Obscuration radius in mm.
        spider = SpiderInfo(
            type=js['Spider']['type']['val'],
            thickness=js['Spider']['thickness']['val'])
        wavelength, efficiency = extTel(js)
        opt_efficiency = OpticalEfficiency(
            wavelength=wavelength, efficiency=efficiency)

        ap_data    = None
        total_area = None
        if spider.type == 'tripod':
            n_apcell = int(epd+4)   # Assuming ap-cell scale to be 1mm/ap-cell.
                                    # Set the aperture pattern size to be 2mm larger than D.
        else:
            raise RuntimeError('Wrong spider type: {}.'.format(spider.type))
        if n_apcell%2 == 1: # n_apcell should be even.
            n_apcell = n_apcell + 1

        ap_data, total_area_mm2 = aperture.calc_aperture(
            n_apcell, epd, r_obscuration, spider.thickness)

        telescope = Telescope(
            epd=epd,
            aperture=ap_data,
            cobs=cobs,
            efl=efl,
            total_area=total_area_mm2*1.e-6,
            spider=spider,
            opt_efficiency=opt_efficiency)
        # total_area is in m^2.

        return telescope


@dataclasses.dataclass(frozen=True)
class Variability():
    """
    Summary: stellar variability class.

    Attributes:
        var     (bool)     : Variability flag.
        Nvar    (int)      : Number of variable objects.
        plate   (list[int]): Plate number for each object.
        star    (list[int]): Star ID number for each object.
        vartype (list[str]): Variable type.
        varfile (list[str]): File name which contains variable template.
        dirname (list[str]): Directory which contains `varfile`.
    """
    var    : bool
    Nvar   : int
    plate  : List[int] = dataclasses.field(default_factory=list)
    star   : List[int] = dataclasses.field(default_factory=list)
    vartype: List[str] = dataclasses.field(default_factory=list)
    varfile: List[str] = dataclasses.field(default_factory=list)
    dirname: List[str] = dataclasses.field(default_factory=list)

    @classmethod
    def from_json(cls,var_json_filename):
        """
        Summary: this function is json i/o for stellar variability
        """
        with open(var_json_filename, "r") as f:
            var_params = json.load(f)
        Nvar=len(var_params)
        plate=[]
        star=[]
        vartype=[]
        varfile=[]
        dirname=[]
        for i in range(0,Nvar):
            var=var_params[str(i)]
            plate.append(var["plate"])
            star.append(var["star"])
            vartype.append(var["vartype"])
            varfile.append(var["varfile"])
            dirname.append(var["dirname"])
        var = Variability(
            var=True, Nvar=Nvar, plate=plate, star=star,
            vartype=vartype, varfile=varfile, dirname=dirname)
        return var


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
                    injlc, b=transitmodel.gentransit_json(
                        t_day,os.path.join(self.dirname[i],self.varfile[i]))
                    sw=True
                    return sw,injlc,b
                else:
                    sys.exit("No valid vartype")

        sw=False
        return sw,None,None


@dataclasses.dataclass(frozen=False)
class Drift:
    """
    Summary: drift class.

    Attributes:
        dft            (bool)  : Drift flag.
        drift_velocity (float) : Drift velocity in ???.
        drift_azimuth  (float) : Drift azimuth angle in ???.
        drift_time     (float) : Drift duration in second.
        drift_length   (float) : ??? length of the drift motion in ???.
        drift_theta    (float) : trajectory?
    """
    dft           : bool
    drift_velocity: float
    drift_azimuth : float
    drift_time    : float = dataclasses.field(init=False)
    drift_length  : float = dataclasses.field(init=False)
    drift_theta   : float = dataclasses.field(init=False)

    @classmethod
    def from_json(self,dft_json_filename):
        """
        Summary: this function is json i/o for drift
        """
        with open(dft_json_filename, "r") as f:
            var_params = json.load(f)
        velocity = var_params["linear"]["drift_velocity"]
        azimuth  = var_params["linear"]["drift_azimuth"]
        dft = Drift(dft=True, drift_velocity=velocity, drift_azimuth=azimuth)
        return dft

    def compute_drift(self,dtace,Nace):
        from jis.pixsim import gentraj

        self.drift_time=dtace*Nace
        self.drift_length=self.drift_velocity*self.drift_time
        self.drift_theta=gentraj.gentraj_drift(Nace,self.drift_length,self.drift_azimuth)
