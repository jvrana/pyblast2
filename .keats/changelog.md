# pyblastbio change log
## 0.5.1
**2019-11-16T11:35:39.742728**
refactoring, reformatting, bug fixes

 - TempFile model ensures temporary files are properly closed and removed
 - Import buf fixes


## 0.5
**2019-10-18T12:42:18.611979**
bug fixed for API breaking changes to networkx 2.4

 - update networkx to networkx 2.4
 - update deps


## 0.4.0a1
**2019-09-30T16:27:15.102015**
span initialization changes




## 0.4.1
**2019-09-30T14:11:10.016367**
minor feature

 - length now in Span.__repr__


## 0.4.0
**2019-09-30T14:01:36.175177**
major revisions to result parsing

 - fully wrapped circular alignments are now supported. The full alignment can be accessed from `x["query"]["start"], x["query"]["raw_end"]`
 - new `result_to_span` function available in the BioBlast class. Function converts a restult to a Span instance (no direction available).
 - new `reindex` method for Span. Method will reindex and adjust the `a`, `b`, and `c` attributes
 - Spans may now wrap around context multiple times. Mapped endpoint accessible by `b` attr and unmapped endpoint accessible by `c` attribute.


## 0.3.8
**2019-09-15T12:08:52.506445**
new features

 - Span now has option to specify full wrapped regions


## 0.3.7
**2019-09-14T08:51:11.408264**
update typing




## 0.3.6
**2019-09-13T07:55:06.099602**
update deps

 - loggable-jdv updated to ^0.1.4 which fixes bug with raising exceptions


## 0.3.5
**2019-09-12T11:14:01.715031**
update deps




## 0.3.4
**2019-09-03T20:54:34.990067**
bug fixes

 - correct edge case bug in Sspan.new(L, L) when span length is L and span is cyclic (previously would return the entire span, now return span(0, 0, L))


## 0.3.3
**2019-09-03T08:27:06.737215**
bug fix

 - fixes Span.spans_origin bug


## 0.3.2
**2019-09-01T19:42:57.784247**
bug fixes

 - fixes initializion bug that would improperly set Span.b if b wrapped around origin


## 0.3.1
**2019-09-01T19:42:09.566888**
poetry bump




## 0.2.17
**2019-08-31T07:39:26.219253**
new features

 - removal of temporary files to avoid error after running blast many times


## 0.2.16
**2019-08-20T14:00:31.145337**
api changes

 - adding records to bioblastfactory returns (keys, records)
 - bioblastfactory now holds config, which is instantiated in the blast when new instance is created
 - blast.__init__ uses config as argument instead of **kwargs
 - BioBlastFactory.add_records returns (keys, records) instead of just records


## 0.2.15
**2019-08-19T11:49:45.194774**
bug fixes

 - fixes .sub index bugs


## 0.2.14
**2019-08-18T08:31:01.427353**
bug fixes

 - fixes bugs with .sub span


## 0.2.13
**2019-08-18T08:13:52.406647**
bug fixes

 - raises error for invalid span over origin


## 0.2.12
**2019-08-17T14:46:07.468458**
bug fixes to span

 - subspan throws appropriate errors for invalid indices


## 0.2.11
**2019-08-16T15:09:55.078414**
non-api changing

 - removed CMD print statement


## 0.2.10
**2019-08-13T20:17:29.764489**
serious bug fixes

 - indices of alignments over spans now fixed (off-by-one error)
 - reversed alignments over spans now properly return their start and ending positions


## 0.2.9
**2019-08-13T19:01:18.755123**
bug fixes




## 0.2.8
**2019-08-08T08:27:05.972192**
bug fix

 - filters out self alignments (alignments with same origin_key and identical end points)


## 0.2.7
**2019-08-07T22:52:21.935282**
new features

 - added BioBlastFactory
 - removed TmpBlast and BlastBase from main import route. These are accessible via `pyblast.blast.TmpBlast` etc.
 - records are now forced to have unique record ids
 - `force_unique_record_id` method provided in `pyblast.utils`
 - SeqRecordDB can only apply one type of transformation once per record key
 - fixed bugs associated with circular subjects and queries
 - Start and ends now return indices adjusted back to origin record if record was pseudocircularized
 - Remove 'alignments' as SeqFeatures are incompatible with circular sequences.


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
