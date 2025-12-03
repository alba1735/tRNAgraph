import sys
import os
import unittest
from unittest.mock import MagicMock, patch
from types import SimpleNamespace

# Add current directory to path
sys.path.append(os.getcwd())

import toolsMap

class TestMapSamples(unittest.TestCase):
    def test_default_bamdir(self):
        args = SimpleNamespace(
            database="db",
            experiment="test_exp",
            samples="samples.txt",
            gtf=None,
            bed=None,
            lazy=False,
            nofrag=False,
            nosizefactors=False,
            bamdir=None,
            threads=1,
            minnontrnasize=20,
            local=False,
            skipcheck=False,
            maxmismatches=None,
            mincoverage=None,
            uniqueonly=False,
            pairs=None,
            paironly=False,
            hubonly=False,
            hub=False,
            maponly=False,
            dumpother=False
        )
        
        mapper = toolsMap.MapSamples(args)
        self.assertEqual(mapper.bamdir, "bam/test_exp")
        
    def test_explicit_bamdir(self):
        args = SimpleNamespace(
            database="db",
            experiment="test_exp",
            samples="samples.txt",
            gtf=None,
            bed=None,
            lazy=False,
            nofrag=False,
            nosizefactors=False,
            bamdir="custom_bam",
            threads=1,
            minnontrnasize=20,
            local=False,
            skipcheck=False,
            maxmismatches=None,
            mincoverage=None,
            uniqueonly=False,
            pairs=None,
            paironly=False,
            hubonly=False,
            hub=False,
            maponly=False,
            dumpother=False
        )
        
        mapper = toolsMap.MapSamples(args)
        self.assertEqual(mapper.bamdir, "custom_bam")

    @patch('os.makedirs')
    @patch('os.path.exists')
    def test_makedirs_called(self, mock_exists, mock_makedirs):
        args = SimpleNamespace(
            database="db",
            experiment="test_exp",
            samples="samples.txt",
            gtf=None,
            bed=None,
            lazy=False,
            nofrag=False,
            nosizefactors=False,
            bamdir=None,
            threads=1,
            minnontrnasize=20,
            local=False,
            skipcheck=False,
            maxmismatches=None,
            mincoverage=None,
            uniqueonly=False,
            pairs=None,
            paironly=False,
            hubonly=False,
            hub=False,
            maponly=False,
            dumpother=False
        )
        
        # Mock exists to return False so makedirs is called
        mock_exists.return_value = False
        
        # Mock other methods called in main to avoid side effects
        with patch('toolsMap.MapSamples.mapsamples'), \
             patch('toolsMap.MapSamples.makefeaturebed'), \
             patch('toolsMap.MapSamples.countfeatures'), \
             patch('toolsMap.MapSamples.counttypes'), \
             patch('toolsMap.MapSamples.gettrnacoverage'), \
             patch('toolsMap.MapSamples.run_deseq2'):
            
            mapper = toolsMap.MapSamples(args)
            mapper.main()
            
            # Check if os.makedirs was called with bam/test_exp
            mock_makedirs.assert_any_call("bam/test_exp")

if __name__ == '__main__':
    unittest.main()
