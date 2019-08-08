# pyblastbio change log
## 0.2.7
**2019-08-07T19:14:27.645353**
new features

 - added BioBlastFactory
 - removed TmpBlast and BlastBase from main import route. These are accessible via `pyblast.blast.TmpBlast` etc.
 - records are now forced to have unique record ids
 - `force_unique_record_id` method provided in `pyblast.utils`
 - SeqRecordDB can only apply one type of transformation once per record key


## 0.2.6
**2019-08-07T16:49:46.556781**
removed unnecessary dep

 - remove 'lobio' as dep


## 0.2.5
**2019-07-17T07:17:59.657293**
minor change

 - better installation description


## 0.2.4
**2019-07-16T15:49:37.922743**
repackaged




## 0.2.3
**2019-07-16T15:48:19.198387**
deps update

 - updated lobio to 0.5.0


## 0.2.2
**2019-07-16T09:18:03.569377**
updated pkg deps

 - lobio bumped to 0.4.0


## 0.2.0
**2019-07-14T12:44:19.114456**
added seq_db sharing for BioBlast and JSONBlast

 - seq_db option for __init__ of BioBlast and JSONBlast


## 0.1.5
**2019-07-14T11:31:30.535249**
update pkg requirements




## 0.1.4
**2019-07-14T11:25:15.007238**
updated lobio requirement




## 0.1.3
**2019-07-14T08:23:43.148726**
update package deps




## 0.1.2
**2019-07-13T14:13:09.896478**
features

 - tests for short blastn
 - new fasta utility functions
