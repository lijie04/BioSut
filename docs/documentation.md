## Installation

Install biosut through PyPi:
```
pip3 install biosut
```

## usage
### Module biosys
```python3
from biosut.biosys import files,path
```

#### files related operations
For example, check the existance of a file/files and check if it is empty.
```python3
f = 'file.txt'
files.check_exist(f, check_empty=True)
```
Get prefix of a file.

```python3
files.get_prefix(f, include_path=True)
```

#### path related operations
For example, make sure a directory exists.
```python3
o = './test/a'
path.sure_exist(o)
```

or
```python3
new.o = path.sure_exist(o)
```


