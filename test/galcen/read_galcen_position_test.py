from jis.galcen.read_galcen_position import load_jscon_random_stars


def test_load_jscon_random_stars():
    data = load_jscon_random_stars()
    assert len(data["ra"]) == 342157


if __name__ == "__main__":
    test_load_jscon_random_stars()