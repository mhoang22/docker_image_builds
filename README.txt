How to build images on the cluster
More detailed note: https://washu.atlassian.net/wiki/spaces/RUD/pages/1705115761/Docker+and+the+RIS+Compute1+Platform
Basics: https://washu.atlassian.net/wiki/spaces/RUD/pages/1865285780/Docker+Tutorial?atlOrigin=eyJpIjoiNDQ4MTAyMGFjN2MyNDIxOGJhOGNlNWZhNjM1OTMyOTgiLCJwIjoiYyJ9

Step 0 (need to do once, in usernode, not storage): Log in to dockerhub on server (then enter user name and password)

LSB_DOCKER_LOGIN_ONLY=1 \
  bsub -G compute-oncology -q oncology-interactive -Is -a 'docker_build' -- . 


Step 1: docker build and push (the bsub command will build and immediately push to dockerhub)

cd /storage1/fs1/mgriffit/Active/griffithlab/gc2596/hmy/docker_image_builds
cd multiqc

bsub -G compute-oncology -q oncology-interactive -Is -a 'docker_build(myhoang04/multiqc)' -- --tag myhoang04/multiqc:1.21 .

~~~~
cd genomics/

bsub -G compute-oncology -q oncology-interactive -Is -a 'docker_build(myhoang04/genomics)' -- --tag myhoang04/genomics:1.0 --tag myhoang04/genomics:latest .