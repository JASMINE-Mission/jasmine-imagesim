import pytest
import pkg_resources


def load_jscon_data():
    image_ = pkg_resources.resource_filename(
        'jasmine-imagesim', 'data/test_jscon_image.npz')
