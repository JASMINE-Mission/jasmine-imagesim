import pytest
import numpy as np

def generate_subgrids(ix,iy,Nside=25):
    """generate subgrid on detector coordinate 

    Args:
        ix (int): pixel x index
        iy (int): pixel y index
        Nside (int, optional): # of subgrids in one side. Defaults to 25. 
        
    Note:
        Nframe = Nside x Nside

    Returns:
        2d array: epsf_subgrids (2, Nframe)
    """
    epsf_subgrids = np.mgrid[0:Nside,0:Nside]
    epsf_subgrids = epsf_subgrids/Nside
    xy = np.array([ix,iy])
    return epsf_subgrids.reshape(2,Nside*Nside) + xy[:,np.newaxis]


def test_generate_subgrids():
    mg = generate_subgrids(1,2)
    ref = np.array([1.96,2.96])
    assert np.all(mg[:,624]-ref == pytest.approx(0.0))
    assert np.all(np.shape(mg)-np.array([2,625])==0)
    
if __name__ == "__main__":
    test_generate_subgrids()
        