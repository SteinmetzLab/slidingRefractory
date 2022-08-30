# slidingRefractory
Code to perform a new test of whether neurons have contaminated refractory periods, with a sliding window


## Python

### Minimal working example

```python
rpMetrics, cont, rp = slidingRP_all(spikeTimes, spikeClusters, params)
```


### Run tests
```commandline
 pytest python/test_*
```

