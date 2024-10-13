from quspin_core._bindings.array import hello


def test_hello():
    assert hello() == "Hello, World!"
