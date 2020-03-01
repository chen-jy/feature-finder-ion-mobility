"""Common utilities for the LC-IMS-MS/MS feature finder and peak picker.
"""

from typing import Any, List, Tuple

import pyopenms as ms


def get_spectrum_points(spec: ms.MSSpectrum) -> List[List[float]]:
    """Extracts the retention times, mass to charges, intensities, and ion mobility values of all
    peaks in a spectrum.

    Keyword arguments:
    spec: the spectrum to extract data points from

    Returns: a list of lists, where each interior list holds RT, m/z, intensity, and IM data (in
    that order) for a single peak in the spectrum.
    """
    point_data = zip(*spec.get_peaks(), spec.getFloatDataArrays()[0])
    return [[spec.getRT(), mz, intensity, im] for mz, intensity, im in point_data]


def get_im_extrema(spectra: List[ms.MSSpectrum]) -> Tuple[float, float]:
    """Finds the smallest and largest IM values in a list of spectra.

    Keyword arguments:
    spectra: the list of spectra with IM data to scan through

    Returns: a tuple of the smallest and largest IM values, in that order.
    """
    smallest_im, largest_im = float('inf'), float('-inf')

    for spec in spectra:
        points = get_spectrum_points(spec)
        for point in points:
            if point[3] < smallest_im:
                smallest_im = point[3]
            if point[3] > largest_im:
                largest_im = point[3]

    return smallest_im, largest_im


def polygon_area(polygon: List[Tuple[float, float]]) -> float:
    """Computes the area of a convex polygon using the shoelace formula."""
    area = 0.0
    for i in range(len(polygon)):
        area += polygon[i][0] * polygon[(i + 1) % len(polygon)][1]
        area -= polygon[i][1] * polygon[(i + 1) % len(polygon)][0]
    return abs(area) / 2.0


def similar_features(feature1: Any, feature2: Any, rt_threshold: float = 5.0, mz_threshold: float = 0.01) -> bool:
    """Checks if the RTs and m/zs of two features are within fixed thresholds of each other."""
    if isinstance(feature1, ms.Feature) and isinstance(feature2, ms.Feature):
        return (abs(feature1.getRT() - feature2.getRT()) < rt_threshold and
                abs(feature1.getMZ() - feature2.getMZ()) < mz_threshold)
    elif isinstance(feature1, list) and isinstance(feature2, list):
        return (abs(feature1[0] - feature2[0]) < rt_threshold and
                abs(feature1[1] - feature2[1]) < mz_threshold)
    elif isinstance(feature1, ms.Feature) and isinstance(feature2, list):
        return (abs(feature1.getRT() - feature2[0]) < rt_threshold and
                abs(feature1.getMZ() - feature2[1]) < mz_threshold)
    elif isinstance(feature1, list) and isinstance(feature2, ms.Feature):
        return similar_features(feature2, feature1, rt_threshold, mz_threshold)
    else:
        return False


def similar_features_im(feature1: List[float], feature2: List[float], rt_threshold: float = 5.0,
                        mz_threshold: float = 0.01, im_threshold: float = 0.031) -> bool:
    """Checks if the RTs, m/zs, and IMs of two features are within fixed thresholds of each other."""
    return similar_features(feature1, feature2) and abs(feature1[2] - feature2[2]) < im_threshold


def has_peaks(exp: ms.MSExperiment) -> bool:
    """Checks if any spectrum in an experiment has peaks."""
    spectra = exp.getSpectra()
    for spec in spectra:
        if spec.size() > 0:
            return True
    return False
