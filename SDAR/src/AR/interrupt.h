#pragma once

#include "AR/information.h"

namespace AR {

    //! binary interrupt status
    enum class InterruptStatus{change=1, merge=2, destroy=3, none=0,
                               forced_merge=4};

    //! Binary interrupt recoreder
    template <class Tparticle>
    struct InterruptBinary {
        AR::BinaryTree<Tparticle>* adr; // binary tree address
        Float time_now;  // current time
        Float time_end;  // finishing time 
        InterruptStatus status; // binary status
        // Step-limit forced merger request.  It is created at an accepted AR
        // step boundary, identifies the pair by id (the request is consumed by
        // the leaf that carries that pair), and records the expected
        // coalescence time for diagnostics.
        bool forced = false;
        long long forced_id[2] = {-1, -1};
        Float theoretical_merger_time = -1.0;

        InterruptBinary(): adr(NULL), time_now(0.0), time_end(0.0), status(InterruptStatus::none),
                           forced(false), theoretical_merger_time(-1.0) {}

        //InterruptBinary(AR::BinaryTree<Tparticle>* _adr, Tparticle* _p, const Float& _time_now, const Float& _time_end, const InterruptStatus& _status): adr(_adr), time_now(_time_now), time_end(_time_end), status(_status) {}

        //InterruptBinary(InterruptBinary<Tparticle>& _bin) {
        //    adr = _bin.adr;
        //    time_now= _bin.time_now;
        //    time_end= _bin.time_end;
        //    status = _bin.status;
        //}

        void clear() { 
            adr=NULL; 
            time_now=0.0; 
            time_end=0.0; 
            status=InterruptStatus::none;
            forced=false;
            forced_id[0]=forced_id[1]=-1;
            theoretical_merger_time=-1.0;
        }

        bool checkParams() {
            ASSERT((status!=InterruptStatus::none&&adr!=NULL)||status==InterruptStatus::none);
            return true;
        }
    };
}
