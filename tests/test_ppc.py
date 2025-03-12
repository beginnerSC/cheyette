import unittest, subprocess

class TestUtils(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        subprocess.check_call(['poetry', 'build'])

    def test_add_sub(self):
        from ppc import add, sub   # cannot import at the top or otherwise the pyd file is locked so poetry can't build in setUpClass
        self.assertEqual(add(1, 2), 3)
        self.assertEqual(sub(1, 2), -1)

if __name__ == '__main__': 
    unittest.main()