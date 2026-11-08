#pragma once

#ifdef BHMERGER

#include <memory>
#include <vector>
#include <fstream>
#include <cstdlib>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <cerrno>
#include <cstdio>
#include <sys/stat.h>
#include <fcntl.h>
#include <iomanip>
#include <sstream>
#include <string>
#include <unistd.h>

//! Thread- and process-safe shared output for black-hole merger events.
class BHMergerOutput {
public:
    struct ParticleState {
        long long id;
        double mass;
        double pos[3];
        double vel[3];
        double spin[3];
    };

private:
    std::string filename_;
    int precision_;
    bool enabled_;
    // Shared ownership supports ARInteraction copies; each stream has exactly
    // one OpenMP writer. Flush only after the worker region has completed.
    using ChangeStreams = std::vector<std::unique_ptr<std::ofstream>>;
    std::shared_ptr<ChangeStreams> change_streams_;


    static bool writeLocked(const std::string& _filename,
                            const std::string& _content,
                            const bool _truncate,
                            const bool _only_empty = false) {
        const int flags = O_WRONLY | O_CREAT | O_APPEND;
        const int fd = ::open(_filename.c_str(), flags, 0666);
        if (fd < 0) return false;

        struct flock lock;
        lock.l_type = F_WRLCK;
        lock.l_whence = SEEK_SET;
        lock.l_start = 0;
        lock.l_len = 0;
        lock.l_pid = 0;

        while (::fcntl(fd, F_SETLKW, &lock) < 0) {
            if (errno != EINTR) {
                ::close(fd);
                return false;
            }
        }

        struct stat info;
        if (::fstat(fd, &info)<0 || (_truncate && ::ftruncate(fd,0)<0)) {
            ::close(fd);
            return false;
        }
        if (_only_empty && info.st_size>0) {
            ::close(fd);
            return true;
        }
        const off_t original_size = _truncate ? 0 : info.st_size;
        const char* data = _content.data();
        std::size_t remaining = _content.size();
        bool success = true;
        while (remaining > 0) {
            const ssize_t written = ::write(fd, data, remaining);
            if (written > 0) {
                data += written;
                remaining -= static_cast<std::size_t>(written);
            }
            else if (written < 0 && errno == EINTR) {
                continue;
            }
            else {
                success = false;
                break;
            }
        }

        if (!success && S_ISREG(info.st_mode)) {
            // Keep the lock while rolling back an incomplete row.
            if (::ftruncate(fd, original_size)<0)
                std::fprintf(stderr, "Cannot roll back partial merger log row in %s\n", _filename.c_str());
        }
        lock.l_type = F_UNLCK;
        while (::fcntl(fd, F_SETLK, &lock) < 0 && errno == EINTR) {}
        if (::close(fd)<0) success=false;
        return success;
    }

    static void writeParticle(std::ostream& _out,
                              const ParticleState& _particle,
                              const int _width) {
        _out << std::setw(_width) << _particle.id
             << std::setw(_width) << _particle.mass;
        for (int k = 0; k < 3; ++k) {
            _out << std::setw(_width) << _particle.pos[k];
        }
        for (int k = 0; k < 3; ++k) {
            _out << std::setw(_width) << _particle.vel[k];
        }
        for (int k = 0; k < 3; ++k) {
            _out << std::setw(_width) << _particle.spin[k];
        }
    }

    static void writeParticleTitle(std::ostream& _out,
                                   const std::string& _prefix,
                                   const int _width) {
        _out << std::setw(_width) << (_prefix + ".id")
             << std::setw(_width) << (_prefix + ".mass")
             << std::setw(_width) << (_prefix + ".pos.x")
             << std::setw(_width) << (_prefix + ".pos.y")
             << std::setw(_width) << (_prefix + ".pos.z")
             << std::setw(_width) << (_prefix + ".vel.x")
             << std::setw(_width) << (_prefix + ".vel.y")
             << std::setw(_width) << (_prefix + ".vel.z")
             << std::setw(_width) << (_prefix + ".spin.x")
             << std::setw(_width) << (_prefix + ".spin.y")
             << std::setw(_width) << (_prefix + ".spin.z");
    }

public:
    BHMergerOutput(): filename_(), precision_(17), enabled_(false) {}

    template <class Tparticle>
    static ParticleState capture(const Tparticle& _particle) {
        ParticleState state;
        state.id = static_cast<long long>(_particle.id);
        state.mass = _particle.mass;
        state.pos[0] = _particle.pos.x;
        state.pos[1] = _particle.pos.y;
        state.pos[2] = _particle.pos.z;
        state.vel[0] = _particle.vel.x;
        state.vel[1] = _particle.vel.y;
        state.vel[2] = _particle.vel.z;
        state.spin[0] = _particle.spin.x;
        state.spin[1] = _particle.spin.y;
        state.spin[2] = _particle.spin.z;
        return state;
    }

    void open(const std::string& _filename,
              const bool _append,
              const bool _initialize_file,
              const int _width,
              const int _precision) {
        filename_ = _filename;
        precision_ = _precision;
        enabled_ = true;

        // Rank 0 alone creates a fresh file and writes its header. The caller
        // performs an MPI barrier before merger integration starts.
        if (_initialize_file) {
            std::ostringstream header;
            header << std::setw(_width) << "time";
            writeParticleTitle(header, "before1", _width);
            writeParticleTitle(header, "before2", _width);
            header << std::setw(_width) << "binary.ecc"
                   << std::setw(_width) << "binary.semi";
            writeParticleTitle(header, "remnant", _width);
            header << '\n';
            enabled_ = writeLocked(filename_, header.str(), !_append, _append);
        }
    }

    // Change rows use explicit separators: round-trip precision can exceed
    // WRITE_WIDTH, in which case fixed-width columns alone would run together.
    void openChange(const std::string& filename, bool append,
                    bool initialize_file, int precision) {
        change_streams_.reset();
        filename_ = filename;
        precision_ = precision;
        enabled_ = true;
        if (initialize_file)
            enabled_ = writeLocked(filename_,
                "time_old time_new id1 id2 m1 m2 a_old e_old a_new e_new "
                "tmerge_old tmerge_new rx ry rz vx vy vz dynamics_dlog_tgw\n", !append, append);
    }

    void openChangeSharded(const std::string& prefix, bool append,
                           int rank, int threads, int precision) {
        filename_ = prefix;
        precision_ = precision;
        enabled_ = true;
        change_streams_ = std::make_shared<ChangeStreams>();
        for (int t=0; t<threads; ++t) {
            const auto name = prefix + "." + std::to_string(rank) + "." + std::to_string(t);
            std::ifstream existing(name, std::ios::binary | std::ios::ate);
            const bool header = !append || !existing || existing.tellg()==std::streampos(0);
            auto stream = std::make_unique<std::ofstream>(name,
                std::ios::out | (append ? std::ios::app : std::ios::trunc));
            if (!*stream) enabled_ = false;
            if (header) *stream << "time_old time_new id1 id2 m1 m2 a_old e_old a_new e_new "
                "tmerge_old tmerge_new rx ry rz vx vy vz dynamics_dlog_tgw\n";
            *stream << std::setprecision(precision_);
            stream->flush();
            if (!*stream) enabled_ = false;
            change_streams_->push_back(std::move(stream));
        }
        if (threads<=0) enabled_=false;
    }

    void flushChanges() const {
        if (!change_streams_) return;
        for (auto& stream: *change_streams_) {
            stream->flush();
            if (!*stream) { std::fprintf(stderr,"Failed to flush change log\n"); std::abort(); }
        }
    }

    template<class Tparticle>
    void writeChange(double time_old, double time_new,
                     const Tparticle& p1, const Tparticle& p2,
                     double a_old, double e_old, double a_new, double e_new,
                     double tmerge_old, double tmerge_new, double dynamics_dlog_tgw=0.0) const {
        if (!enabled_) return;
        if (change_streams_) {
            int thread=0;
#ifdef _OPENMP
            thread=omp_get_thread_num();
            if (omp_get_level()>1) std::abort(); // thread IDs must be unique writers
#endif
            if (static_cast<std::size_t>(thread)>=change_streams_->size()) std::abort();
            auto& line=*(*change_streams_)[thread];
            line << time_old << ' ' << time_new << ' ' << p1.id << ' ' << p2.id
                 << ' ' << p1.mass << ' ' << p2.mass
                 << ' ' << a_old << ' ' << e_old << ' ' << a_new << ' ' << e_new
                 << ' ' << tmerge_old << ' ' << tmerge_new;
            for (int k=0;k<3;++k) line << ' ' << p1.pos[k]-p2.pos[k];
            for (int k=0;k<3;++k) line << ' ' << p1.vel[k]-p2.vel[k];
            line << ' ' << dynamics_dlog_tgw << '\n';
            if (!line) { std::fprintf(stderr,"Failed to write change log\n"); std::abort(); }
            return;
        }
#pragma omp critical(bhmerger_output)
        {
            std::ostringstream line;
            line << std::setprecision(precision_)
                 << time_old << ' ' << time_new
                 << ' ' << p1.id << ' ' << p2.id
                 << ' ' << p1.mass << ' ' << p2.mass
                 << ' ' << a_old << ' ' << e_old
                 << ' ' << a_new << ' ' << e_new
                 << ' ' << tmerge_old << ' ' << tmerge_new;
            for (int k=0; k<3; ++k) line << ' ' << p1.pos[k]-p2.pos[k];
            for (int k=0; k<3; ++k) line << ' ' << p1.vel[k]-p2.vel[k];
            line << ' ' << dynamics_dlog_tgw << '\n';
            if (!writeLocked(filename_, line.str(), false))
                std::fprintf(stderr, "Failed to write merger event log %s\n", filename_.c_str());
        }
    }

    bool isOpen() const {
        return enabled_;
    }

    void close() {
        flushChanges();
        change_streams_.reset();
        enabled_ = false;
        filename_.clear();
    }

    void write(const double _time,
               const ParticleState& _before1,
               const ParticleState& _before2,
               const double _ecc,
               const double _semi,
               const ParticleState& _remnant,
               const int _width) {
        if (!enabled_) return;

#pragma omp critical(bhmerger_output)
        {
            std::ostringstream line;
            line << std::setprecision(precision_)
                 << std::setw(_width) << _time;
            writeParticle(line, _before1, _width);
            writeParticle(line, _before2, _width);
            line << std::setw(_width) << _ecc
                 << std::setw(_width) << _semi;
            writeParticle(line, _remnant, _width);
            line << '\n';
            if (!writeLocked(filename_, line.str(), false))
                std::fprintf(stderr, "Failed to write merger event log %s\n", filename_.c_str());
        }
    }
};

#endif
