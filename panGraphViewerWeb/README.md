## Welcome to the Web browser-based panGraphViewer

There are 2 ways to install: using ``miniconda3`` or using ``docker``

### 1. Using miniconda
You may use ``pip`` in ``miniconda3`` to install the packages

On ``Linux`` or ``macOS``:
```
pip install -r requirements.txt 
```
To start the server:
```
./run.sh
```

---
On Windows:
```
pip install -r requirements_windows.txt
```
One more packages on ``Window`` platform is
```
conda install m2-base
```
To start the server:
```
run.bat
```

---

And then, you can type ``http://localhost:8000`` or ``http://<IPaddress>:8000`` in your browser to use ``panGraphViewer``

The default demo login info:
```
Account: demo
password: demodemo
```

You can also manage the accounts by using admin account at admin address:
```
Address: http://localhost:8000/admin
Account: admin
password: abcd1234
```

---


### 2. Using docker

If you have docker permission, you can use docker to simplify the deployment.
Go to the directory ``docker``, build the docker image and create the container.

On ``Linux``  or ``macOS``:
```
./build_docker.sh
./run_docker.sh
```

On Windows:
```
build_docker.bat
run_docker.bat
```

When running run_docker.*, you will see the published port number at the end of execution. Then you can use your browser to access the address ```localhost:<PortNumber>``` or  ```<IPaddress>:<PortNumber>```.
E.g. if the published port number is 32793, then you can type in the address ``http://localhost:32793``


The default demo login info:
```
Account: demo
password: demodemo
```

You can also manage the accounts by using admin account at admin address:
```
Address: http://<IPaddress>:<PortNumber>/admin
Account: admin
password: abcd1234
```

---
Enjoy the use of panGraphViewer!
