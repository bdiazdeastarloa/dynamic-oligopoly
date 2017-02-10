% Run 'main_dynolig.m' as a parallel batch job

Cluster=parcluster('LionX');
Cluster.ResourceTemplate='-l nodes=24 -l walltime=20:00:00 -l pmem=5gb -q lionxf-econ';
job = batch(Cluster,'main','Matlabpool',23,'CaptureDiary',true);
