"""Pytest for utilities functions"""
import numpy as np
import numpy.testing as tst
import utilities.attitude as att

angle = (0 - 2*np.pi) * np.random.random_sample() + 0
x_axis = np.array([1.0,0,0])
y_axis = np.array([0,1.0,0])
z_axis = np.array([0,0,1.0])

pos_90 = np.pi/2
neg_90 = -np.pi/2

rtol = 1e-6
atol = 0
    
def test_ROT1_SO3():
    """Make sure ROT1(angle) \in SO(3)"""

    tst.assert_allclose(np.dot(np.transpose(att.ROT1(angle)),att.ROT1(angle)), np.identity(3))
    tst.assert_allclose(np.linalg.det(att.ROT1(angle)),1.0)

def test_ROT2_SO3():
    """Make sure ROT2(angle) \in SO(3)"""

    tst.assert_allclose(np.dot(np.transpose(att.ROT2(angle)),att.ROT2(angle)), np.identity(3))
    tst.assert_allclose(np.linalg.det(att.ROT2(angle)),1.0)

def test_ROT3_SO3():
    """Make sure ROT1(angle) \in SO(3)"""

    tst.assert_allclose(np.dot(np.transpose(att.ROT3(angle)),att.ROT3(angle)), np.identity(3))
    tst.assert_allclose(np.linalg.det(att.ROT3(angle)),1.0)

def test_ROT1_90():
    """Make sure a 90 degree rotation rotates a column vector correctly"""

    # np.pi/2 rotations for each axis
    tst.assert_array_almost_equal(np.dot(x_axis,att.ROT1(pos_90)),x_axis)
    tst.assert_array_almost_equal(np.dot(y_axis,att.ROT1(pos_90)),z_axis)
    tst.assert_array_almost_equal(np.dot(z_axis,att.ROT1(pos_90)),-y_axis)

    # negative 90 rotations
    tst.assert_array_almost_equal(np.dot(x_axis,att.ROT1(neg_90)),x_axis)
    tst.assert_array_almost_equal(np.dot(y_axis,att.ROT1(neg_90)),-z_axis)
    tst.assert_array_almost_equal(np.dot(z_axis,att.ROT1(neg_90)),y_axis)

def test_ROT2_90():
    """Make sure a 90 degree rotation rotates a column vector correctly"""

    # np.pi/2 rotations for each axis
    tst.assert_array_almost_equal(np.dot(x_axis,att.ROT2(pos_90)),-z_axis)
    tst.assert_array_almost_equal(np.dot(y_axis,att.ROT2(pos_90)),y_axis)
    tst.assert_array_almost_equal(np.dot(z_axis,att.ROT2(pos_90)),x_axis)

    # negative 90 rotations
    tst.assert_array_almost_equal(np.dot(x_axis,att.ROT2(neg_90)),z_axis)
    tst.assert_array_almost_equal(np.dot(y_axis,att.ROT2(neg_90)),y_axis)
    tst.assert_array_almost_equal(np.dot(z_axis,att.ROT2(neg_90)),-x_axis)

def test_ROT3_90():
    """Make sure a 90 degree rotation rotates a column vector correctly"""

    # np.pi/2 rotations for each axis
    tst.assert_array_almost_equal(np.dot(x_axis,att.ROT3(pos_90)),y_axis)
    tst.assert_array_almost_equal(np.dot(y_axis,att.ROT3(pos_90)),-x_axis)
    tst.assert_array_almost_equal(np.dot(z_axis,att.ROT3(pos_90)),z_axis)

    # negative 90 rotations
    tst.assert_array_almost_equal(np.dot(x_axis,att.ROT3(neg_90)),-y_axis)
    tst.assert_array_almost_equal(np.dot(y_axis,att.ROT3(neg_90)),x_axis)
    tst.assert_array_almost_equal(np.dot(z_axis,att.ROT3(neg_90)),z_axis)