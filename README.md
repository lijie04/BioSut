# biosut
**
biosut ("biology suite tool") is a packge containing some biology-related bioinformatics operations on sequences, files, directories. 
**

## Installation
Install biosut through PyPi:
```
pip3 install biosut
```

## usage
### Module biosys
```
from biosut.biosys import files,path
```
#### files related operations
For example, check the existance of a file/files and check if it's empty.
```
f = 'file.txt'
files.check_exist(f, check_empty=True)
```
Get prefix of a file.
```
files.get_prefix(f, include_path=True)
```

#### path related operations
For example, make sure a directory exists.
```
o = './test/a'
path.sure_exist(o)
```
or
```
new.o = path.sure_exist(o)
```



