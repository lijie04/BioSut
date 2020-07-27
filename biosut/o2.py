
from test_object import test

t = test(10, 8)
t.run()
print(t.a)

from biosys import path

cmd = 'samtools view -H %s|head -n 1'
out, err = path.exe_proc(cmd)
print('out', out)
print('err', err)


