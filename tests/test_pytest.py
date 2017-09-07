import pytest


@pytest.fixture
def foo():
    return 5
def test():

    print(foo())
    # https: // docs.pytest.org / en / latest / builtin.html
    # def test1():
    #     with pytest.raises(ValueError) as e:
    #         dosomething()
    #     pytest.approx
    #     