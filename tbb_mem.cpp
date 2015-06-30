#include "tbb_mem.h"

#include "tbb/task_scheduler_init.h"

extern "C" void mem_process_seqs_tbb(ktp_aux_t* aux, int n_threads)
{
	tbb::task_scheduler_init tbb_init(n_threads);

	bwa::ReadBSeqFilter read_bseq(*aux);
	bwa::MemParallelFilter mem(*aux);
	bwa::WriteBSeqFilter write_bseq(*aux);

	tbb::parallel_pipeline(n_threads*2,
		tbb::make_filter<void,bwa::BSeq*>(tbb::filter::serial_in_order,read_bseq)&
		tbb::make_filter<bwa::BSeq*,bwa::BSeq*>(tbb::filter::parallel,mem)&
#ifndef TBB_OUT_OF_ORDER_WRITE
		tbb::make_filter<bwa::BSeq*,void>(tbb::filter::serial_in_order,write_bseq));
#else
		tbb::make_filter<bwa::BSeq*,void>(tbb::filter::serial_out_of_order,write_bseq));
#endif
}
