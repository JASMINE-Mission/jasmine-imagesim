import pytest
from jis.galcen.read_galcen_position import load_jscon_data

def test_load_jscon_data():
    gal_l, gal_b, hw, hwtarget = load_jscon_data()
    assert len(gal_l) == 391359


if __name__ == "__main__":
    load_jscon_data()