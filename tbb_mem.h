#pragma once

#include <stdint.h>
#include <vector>
#include <iostream>
#include <memory>

#include "tbb/pipeline.h"
#include "tbb/tick_count.h"

#ifdef USE_MALLOC_WRAPPERS
#  include "malloc_wrap.h"
#endif
#include "bwamem.h"
#include "bwa.h"
#include "kseq.h"
#include "utils.h"
extern "C"
{
	KSEQ_DECLARE(gzFile)
}

//types are copy/past from fastmap.c
typedef struct {
	kseq_t *ks, *ks2;
	mem_opt_t *opt;
	mem_pestat_t *pes0;
	int64_t n_processed;
	int copy_comment, actual_chunk_size;
	bwaidx_t *idx;
} ktp_aux_t;

typedef struct {
	ktp_aux_t *aux;
	int n_seqs;
	bseq1_t *seqs;
} ktp_data_t;

extern "C" ktp_data_t* process_0(ktp_aux_t* aux);
extern "C" void process_1(ktp_aux_t* aux, ktp_data_t* data);
extern "C" void process_2(ktp_aux_t* aux, ktp_data_t* data);

namespace bwa
{
	class BSeq
	{
	public:
		BSeq(ktp_data_t* data, uint64_t n_processed): data_(data), n_processed_(n_processed)
		{
		}

		~BSeq()
		{
		}

		ktp_data_t* ktp_data()
		{
			return data_;
		}

		uint64_t n_processed()
		{
			return n_processed_;
		}
	
	private:
		ktp_data_t* data_;
		const size_t n_processed_;
	};

	class ReadBSeqFilter
	{
	public:
		ReadBSeqFilter(ktp_aux_t& aux): aux_(aux), n_processed_(0)
		{
			if (bwa_verbose >= 4) std::clog<<"[D::bwa::ReadBSeqFilter::ReadBSeqFilter()]\n";
		}

		~ReadBSeqFilter()
		{
			if (bwa_verbose >= 4) std::clog<<"[D::bwa::ReadBSeqFilter::~ReadBSeqFilter()]\n";
		}

		BSeq* operator()(tbb::flow_control& fc) const
		{
			tbb::tick_count t;
			
			if(bwa_verbose>=3)
			{
				t=tbb::tick_count::now();
				if (bwa_verbose >= 4) std::clog<<"[D::bwa::ReadBSeqFilter::operator()]\n";
			}

			ktp_data_t* data=0;
			if(!(data=process_0(&aux_)))
			{
				fc.stop();
				return 0;
			}
			size_t n_processed=n_processed_;
			n_processed_+=data->n_seqs;;
			
			if(bwa_verbose>=3)	std::clog<<"[M::bwa::ReadBSeqFilter::operator()] "<<(tbb::tick_count::now()-t).seconds()<<" sec\n";
			
			return new BSeq(data,n_processed);
		}
	
	private:
		ktp_aux_t& aux_;
		mutable size_t n_processed_;
	};

	class WriteBSeqFilter
	{
	public:
		WriteBSeqFilter(ktp_aux_t& aux): aux_(aux)
		{
			if (bwa_verbose >= 4) std::clog<<"[D::bwa::WriteBSeqFilter::WriteBSeqFilter()]\n";
		}

		~WriteBSeqFilter()
		{
			if (bwa_verbose >= 4) std::clog<<"[D::bwa::WriteBSeqFilter::~WriteBSeqFilter()]\n";
		}

		void operator()(BSeq* data) const
		{
			tbb::tick_count t;
			
			if(bwa_verbose>=3)
			{
				t=tbb::tick_count::now();
				if (bwa_verbose >= 4) std::clog<<"[D::bwa::WriteBSeqFilter::operator()]\n";
			}

			process_2(&aux_,data->ktp_data());

			delete data;

			if(bwa_verbose>=3)	std::clog<<"[M::bwa::WriteBSeqFilter::operator()] "<<(tbb::tick_count::now()-t).seconds()<<" sec\n";
		}
	
	private:
		ktp_aux_t& aux_;
	};


	class MemParallelFilter
	{
	public:
		MemParallelFilter(ktp_aux_t& aux): aux_(aux)
		{
			if (bwa_verbose >= 4) std::clog<<"[D::bwa::MemParallelFilter::MemParallelFilter()]\n";
		}

		~MemParallelFilter()
		{
			if (bwa_verbose >= 4) std::clog<<"[D::bwa::MemParallelFilter::~MemParallelFilter()]\n";
		}

		BSeq* operator()(BSeq* data) const
		{
			tbb::tick_count t;

			if(bwa_verbose>=3)
			{
				t=tbb::tick_count::now();
				if (bwa_verbose >= 4) std::clog<<"[D::bwa::MemParallelFilter::operator()]\n";
			}

			static ktp_aux_t aux=aux_;
			aux.n_processed=data->n_processed();
			process_1(&aux,data->ktp_data());

			if(bwa_verbose>=3)	std::clog<<"[M::bwa::MemParallelFilter::operator()] "<<(tbb::tick_count::now()-t).seconds()<<" sec\n";

			return data;
		}

	private:
		ktp_aux_t& aux_;
	};
}

