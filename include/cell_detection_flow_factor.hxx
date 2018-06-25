#ifndef LP_MP_CELL_DETECTION_FLOW_FACTOR_HXX
#define LP_MP_CELL_DETECTION_FLOW_FACTOR_HXX

// model the cell detections with a flow network. Divisions are modelled through duplicate duplicated edges.
/*
       /----o
  o----
       \----o

is transformed to

       /----o
  o----

  s---------o         ---t
   \-----------------/

   where s and t are additional nodes. An equality constraint (enforced by other factors) ensures that flow values are equal: o-o = s-o

*/

#include "mcf_ssp.hxx"
#include "cell_tracking_constructor.hxx"

namespace LP_MP {

class cell_detection_flow_factor {
public:
    cell_detection_flow_factor(const std::size_t no_cells, const std::size_t no_cell_transitions, const std::size_t no_cell_divisions)
    {
        const std::size_t no_mcf_nodes = 2 + 2*no_cells;
        // additional edges due to (dis)appearances, cell divisions, source to sink
        const std::size_t no_mcf_edges = 1 + 3*no_cells + no_cell_transitions + 2*no_cell_divisions;
        mcf_ = MCF::SSP<int,REAL>(no_mcf_nodes, no_mcf_edges); 

        mcf_.add_node_excess(0, no_cells);
        mcf_.add_node_excess(no_mcf_nodes-1, -int(no_cells));
        mcf_.add_edge(0, no_mcf_nodes-1, 0, no_cells, -1e-8);
    }

    // returns appearance and disappearance edge
    std::array<std::size_t,3> add_detection_hypothesis(const std::size_t cell_detection_counter, const double detection_cost, const double appearance_cost, const double disappearance_cost)
    {
        assert(detection_cost == 0.0 && appearance_cost == 0.0 && disappearance_cost == 0.0);

        const auto appearance_edge = mcf_.add_edge(0, 2*cell_detection_counter + 1, 0, 1, appearance_cost);
        const auto detection_edge = mcf_.add_edge(2*cell_detection_counter + 1, 2*cell_detection_counter + 2, 0, 1, detection_cost);
        const auto disappearance_edge = mcf_.add_edge(2*cell_detection_counter + 2, mcf_.no_nodes()-1, 0, 1, disappearance_cost);

        assert(cell_detection_counter < mcf_.no_nodes()-2);
        return {appearance_edge, detection_edge, disappearance_edge};
    }

    std::size_t add_cell_transition(const std::size_t cell_detection_counter_out, const std::size_t cell_detection_counter_in, const double cost)
    {
        assert(cost == 0.0);
        return mcf_.add_edge(2*cell_detection_counter_out + 2, 2*cell_detection_counter_in + 1, 0, 1, cost);
    }

    std::array<std::size_t,2> add_cell_division(const std::size_t cell_detection_counter_out, const std::size_t cell_detection_counter_in_1, const std::size_t cell_detection_counter_in_2, const double cost)
    {
        assert(cost == 0.0);
        const std::size_t edge_1 = mcf_.add_edge(2*cell_detection_counter_out + 2, 2*cell_detection_counter_in_1 + 1, 0, 1, 0.5*cost);
        const std::size_t edge_2 = mcf_.add_edge(0, 2*cell_detection_counter_in_2 + 1, 0, 1, 0.5*cost); // copied edge
        return {edge_1,edge_2};
    } 

    REAL LowerBound() const
    {
        //std::cout << "mcf lower bound = " << mcf_.solve() << "\n";
        //std::cout << "print flow\n";
        //mcf_.print_flow();
        return mcf_.solve();
    }

    REAL EvaluatePrimal() const
    {
        return std::numeric_limits<REAL>::infinity();
    }

    void MaximizePotential()
    {
        mcf_.solve();
    }
    void MaximizePotentialAndComputePrimal()
    {
        mcf_.solve();
    }

    void init_primal() {}
    auto export_variables() { return std::tie(); }
    template<class ARCHIVE> void serialize_primal(ARCHIVE& ar) { ar(); }
    template<class ARCHIVE> void serialize_dual(ARCHIVE& ar) { ar(); }


    template<typename SOLVER>
    void construct_constraints(SOLVER& s) const
    {}

    template<typename SOLVER>
    void convert_primal(SOLVER& s)
    {}


    mutable MCF::SSP<int,REAL> mcf_; // make private and make cell_detection_flow_message a friend?
private:
}; 

// left is detection factor, right is flow factor
class cell_detection_flow_message {
public:
    template<typename INCOMING_EDGE_ITERATOR, typename OUTGOING_EDGES_ITERATOR>
    cell_detection_flow_message(const std::size_t _detection_edge, INCOMING_EDGE_ITERATOR incoming_begin, INCOMING_EDGE_ITERATOR incoming_end, OUTGOING_EDGES_ITERATOR outgoing_begin, OUTGOING_EDGES_ITERATOR outgoing_end)
    : detection_edge(_detection_edge),
    incoming_edge_indices(incoming_begin, incoming_end),
    outgoing_edge_indices(outgoing_begin, outgoing_end)
    {
#ifndef NDEBUG
        std::vector<std::size_t> edges;
        std::copy(incoming_begin, incoming_end, std::back_inserter(edges));
        for(auto it=outgoing_begin; it!=outgoing_end; ++it) {
            edges.push_back((*it)[0]);
            if((*it)[0] != (*it)[1]) {
                edges.push_back((*it)[1]);
            }
        }
        assert(HasUniqueValues(edges));
#endif
    }

    template<typename RIGHT_FACTOR>
    void RepamRight(RIGHT_FACTOR& r, const REAL msg, const INDEX msg_dim)
    {
        assert(msg_dim < 1 + incoming_edge_indices.size() + 2*outgoing_edge_indices.size());

        if(msg_dim == 0) {
            r.mcf_.update_cost(detection_edge, msg);
        } else if(msg_dim < 1+incoming_edge_indices.size()) {
            r.mcf_.update_cost(incoming_edge_indices[msg_dim-1], msg);
        } else {
            const auto edge = (msg_dim - 1 - incoming_edge_indices.size()) / 2;
            const auto div = (msg_dim - 1 - incoming_edge_indices.size()) % 2;
            r.mcf_.update_cost(outgoing_edge_indices[edge][div], msg);
        } 
    }

    template<typename LEFT_FACTOR>
    void RepamLeft(LEFT_FACTOR& l, const REAL msg, const INDEX msg_dim)
    {
        assert(msg_dim < 1 + incoming_edge_indices.size() + 2*outgoing_edge_indices.size());
        assert(incoming_edge_indices.size() == l.incoming.size() && outgoing_edge_indices.size() == l.outgoing.size());

        if(msg_dim == 0) {
            l.detection += msg;
        } else if(msg_dim < 1+incoming_edge_indices.size()) {
            l.update_incoming(msg_dim-1, msg); 
        } else {
            const auto edge = (msg_dim - 1 - incoming_edge_indices.size()) / 2;
            l.update_outgoing(edge, msg); 
        }
    }

    template<typename LEFT_FACTOR, typename MSG>
    void RepamLeft(LEFT_FACTOR& l, const MSG& msg)
    {
        assert(false); exit(1);
        //std::cout << "msg left = ";
        //for(std::size_t i=0; i<msg.size(); ++i) { std::cout << msg[i] << ", "; }
        //std::cout << "\n";

        l.detection += msg[0];

        assert(msg.size() == 1 + l.incoming.size() + l.outgoing.size());
        assert(incoming_edge_indices.size() == l.incoming.size() && outgoing_edge_indices.size() == l.outgoing.size());

        for(std::size_t i=0; i<incoming_edge_indices.size(); ++i) {
            l.update_incoming(i, msg[1+i]);
        } 

        for(std::size_t i=0; i<outgoing_edge_indices.size(); ++i) {
            l.update_outgoing(i, msg[1+incoming_edge_indices.size() + i]);
        } 
    }

    template<typename RIGHT_FACTOR, typename MSG>
    void RepamRight(RIGHT_FACTOR& r, const MSG& msg)
    {
        assert(false); exit(1);
        assert(msg.size() == 1 + incoming_edge_indices.size() + outgoing_edge_indices.size());

        for(INDEX i=0; i<msg.size(); ++i) {
            if(!std::isfinite(msg[i])) { 
                //throw std::runtime_error("msg out of bound"); 
            }
        }
#ifndef NDEBUG
        const auto left_node = r.mcf_.tail(detection_edge);
        const auto right_node = r.mcf_.head(detection_edge);
        assert(msg.size() == 1 + r.mcf_.no_outgoing_arcs(left_node) + r.mcf_.no_outgoing_arcs(right_node) - 2);
        for(std::size_t i=0; i<incoming_edge_indices.size(); ++i) { assert(r.mcf_.head(incoming_edge_indices[i]) == left_node); }
        for(std::size_t i=0; i<incoming_edge_indices.size(); ++i) { assert(r.mcf_.tail(incoming_edge_indices[i]) < left_node); }
        for(std::size_t i=0; i<outgoing_edge_indices.size(); ++i) { assert(r.mcf_.tail(outgoing_edge_indices[i][0]) == right_node); }
        for(std::size_t i=0; i<outgoing_edge_indices.size(); ++i) { assert(r.mcf_.head(outgoing_edge_indices[i][0]) > right_node); }
        for(std::size_t i=0; i<outgoing_edge_indices.size(); ++i) { assert(r.mcf_.head(outgoing_edge_indices[i][1]) > right_node); }

        assert(std::isfinite(r.mcf_.cost(detection_edge)));
        for(std::size_t i=0; i<incoming_edge_indices.size(); ++i) { assert(std::isfinite(r.mcf_.cost(incoming_edge_indices[i]))); }
        for(std::size_t i=0; i<outgoing_edge_indices.size(); ++i) { assert(std::isfinite(r.mcf_.cost(outgoing_edge_indices[i][0]))); }
        for(std::size_t i=0; i<outgoing_edge_indices.size(); ++i) { assert(std::isfinite(r.mcf_.cost(outgoing_edge_indices[i][1]))); }
#endif
        //std::cout << "msg right = ";
        //for(std::size_t i=0; i<msg.size(); ++i) { std::cout << msg[i] << ", "; }
        //std::cout << "\n";

        r.mcf_.update_cost(detection_edge, msg[0]);

        // first messages are incoming
        for(std::size_t i=0; i<incoming_edge_indices.size(); ++i) {
            const auto edge_idx = incoming_edge_indices[i];
            assert(r.mcf_.lower_bound(edge_idx) == 0 && r.mcf_.upper_bound(edge_idx) == 1);
            r.mcf_.update_cost(edge_idx, msg[1+i]);
        } 

        // then come outgoing edges.
        for(std::size_t i=0; i<outgoing_edge_indices.size(); ++i) {
            const auto out_edge_idx_1 = outgoing_edge_indices[i][0]; 
            const auto out_edge_idx_2 = outgoing_edge_indices[i][1]; 
            assert(r.mcf_.lower_bound(out_edge_idx_1) == 0 && r.mcf_.upper_bound(out_edge_idx_1) == 1);
            assert(r.mcf_.lower_bound(out_edge_idx_2) == 0 && r.mcf_.upper_bound(out_edge_idx_2) == 1);
            if(out_edge_idx_1 == out_edge_idx_2) {
                r.mcf_.update_cost(out_edge_idx_1, msg[1 + i + incoming_edge_indices.size()]); // transition
            } else { // division
                assert(false); // turn on again for general problems with division
                r.mcf_.update_cost(out_edge_idx_1, 0.5*msg[1 + i + incoming_edge_indices.size()]);
                r.mcf_.update_cost(out_edge_idx_2, 0.5*msg[1 + i + incoming_edge_indices.size()]);
            }
        }

        assert(std::isfinite(r.mcf_.cost(detection_edge)));
        for(std::size_t i=0; i<incoming_edge_indices.size(); ++i) { assert(std::isfinite(r.mcf_.cost(incoming_edge_indices[i]))); }
        for(std::size_t i=0; i<outgoing_edge_indices.size(); ++i) { assert(std::isfinite(r.mcf_.cost(outgoing_edge_indices[i][0]))); }
        for(std::size_t i=0; i<outgoing_edge_indices.size(); ++i) { assert(std::isfinite(r.mcf_.cost(outgoing_edge_indices[i][1]))); }
    }

    template<typename LEFT_FACTOR, typename MSG>
    void send_message_to_right(LEFT_FACTOR& l, MSG& msg, const REAL omega)
    {
        /*
        vector<REAL> m(1 + l.incoming.size() + l.outgoing.size());
        m[0] = l.detection;
        for(std::size_t i=0; i<l.incoming.size(); ++i) { m[1 + i] = l.incoming[i]; }
        for(std::size_t i=0; i<l.outgoing.size(); ++i) { m[1 + l.incoming.size() + i] = l.outgoing[i]; }
        for(auto& x : m) { x *= omega; }
        msg -= m;
        return;
        */

        const auto min_incoming_cost = l.min_incoming();
        const auto min_outgoing_cost = l.min_outgoing();
        // possibly think about second min as well
        const auto min_cost = l.min_detection_cost();
        msg[0] -= omega*min_cost;
        for(std::size_t i=0; i<l.incoming.size(); ++i) {
            msg[1 + i] -= omega*(l.incoming[i] - min_incoming_cost);
        }
        for(std::size_t i=0; i<l.outgoing.size(); ++i) {
            const auto msg_val = omega*(l.outgoing[i] - min_outgoing_cost);
            msg[1 + l.incoming.size() + 2*i] -= 0.5*msg_val;
            msg[1 + l.incoming.size() + 2*i+1] -= 0.5*msg_val;
        }
        /*
        vector<REAL> m(1 + l.incoming.size() + l.outgoing.size());
        m[0] = l.detection - min_cost;
        for(std::size_t i=0; i<l.incoming.size(); ++i) { m[1 + i] = l.incoming[i] - min_incoming_cost; }
        for(std::size_t i=0; i<l.outgoing.size(); ++i) { m[1 + l.incoming.size() + i] = l.outgoing[i] - min_outgoing_cost; }
        //std::cout << "msg = ";
        //for(auto& x : m) { std::cout << x << ", "; }
        //std::cout << "\n";

        for(auto& x : m) { x *= omega; }
        msg -= m;
        */
    }

    template<typename RIGHT_FACTOR, typename MSG>
    void send_message_to_left(RIGHT_FACTOR& r, MSG& msg, const REAL omega)
    {
        assert(false);
    }

    template<typename RIGHT_FACTOR, typename MSG_ITERATOR>
    static void SendMessagesToLeft(RIGHT_FACTOR& r, MSG_ITERATOR msg_begin, MSG_ITERATOR msg_end, const REAL omega)
    {
        std::vector<REAL> cover_number(r.mcf_.no_edges(), 0.0); // no of times the message is covered
        for(auto msg_it=msg_begin; msg_it!=msg_end; ++msg_it) {
            auto& msg = (*msg_it).GetMessageOp();
            assert( cover_number[ msg.detection_edge ] == 0);
            cover_number[ msg.detection_edge ]++;
            for(std::size_t i=0; i<msg.incoming_edge_indices.size(); ++i) {
                cover_number[ msg.incoming_edge_indices[i] ]++;
            }
            for(std::size_t i=0; i<msg.outgoing_edge_indices.size(); ++i) {
                if(msg.outgoing_edge_indices[i][0] == msg.outgoing_edge_indices[i][1]) {
                    cover_number[ msg.outgoing_edge_indices[i][0] ]++;
                } else {
                    cover_number[ msg.outgoing_edge_indices[i][0] ]++;
                    cover_number[ msg.outgoing_edge_indices[i][1] ]++;
                }
            } 
        }
        assert(*std::max_element(cover_number.begin(), cover_number.end()) <= 2);

        for(std::size_t iter=0; iter<1; ++iter) {
            r.mcf_.solve();
            //r.mcf_.distribute_potentials();
            std::vector<REAL> reduced_costs(r.mcf_.no_edges());
            for(std::size_t e=0; e<r.mcf_.no_edges(); ++e) {
                reduced_costs[e] = r.mcf_.reduced_cost(e); 
                /*
                if(r.mcf_.flow(e) == r.mcf_.lower_bound(e)) {
                    assert(r.mcf_.reduced_cost(e) >= -eps);
                } else if(r.mcf_.flow(e) == r.mcf_.upper_bound(e)) {
                    assert(r.mcf_.reduced_cost(e) <= eps);
                } else {
                    assert(std::abs(r.mcf_.reduced_cost(e)) <= eps);
                }
                */
            }

            for(auto msg_it=msg_begin; msg_it!=msg_end; ++msg_it) {
                auto& msg = (*msg_it).GetMessageOp();
                //std::cout << &*msg_it << "," << &msg << "\n";

                assert(r.mcf_.flow(msg.detection_edge) == 0 || r.mcf_.flow(msg.detection_edge) == 1);
                for(std::size_t i=0; i<msg.incoming_edge_indices.size(); ++i) { assert(r.mcf_.flow(msg.incoming_edge_indices[i]) == 0 || r.mcf_.flow(msg.incoming_edge_indices[i]) == 1); }
                for(std::size_t i=0; i<msg.outgoing_edge_indices.size(); ++i) { assert(r.mcf_.flow(msg.outgoing_edge_indices[i][0]) == 0 || r.mcf_.flow(msg.outgoing_edge_indices[i][0]) == 1); }
                for(std::size_t i=0; i<msg.outgoing_edge_indices.size(); ++i) { assert(r.mcf_.flow(msg.outgoing_edge_indices[i][1]) == 0 || r.mcf_.flow(msg.outgoing_edge_indices[i][1]) == 1); }

                (*msg_it)[0] -= omega * 1.0/cover_number[msg.detection_edge] * reduced_costs[msg.detection_edge];
                for(std::size_t i=0; i<msg.incoming_edge_indices.size(); ++i) {
                    (*msg_it)[1+i] -= omega * 1.0/cover_number[msg.incoming_edge_indices[i]] * reduced_costs[msg.incoming_edge_indices[i]];
                }
                for(std::size_t i=0; i<msg.outgoing_edge_indices.size(); ++i) {
                    if(msg.outgoing_edge_indices[i][0] == msg.outgoing_edge_indices[i][1]) {
                        (*msg_it)[1 + msg.incoming_edge_indices.size() + 2*i] -= omega * 1.0/cover_number[msg.outgoing_edge_indices[i][0]] * reduced_costs[msg.outgoing_edge_indices[i][0]];
                    } else {
                        assert(cover_number[msg.outgoing_edge_indices[i][0]] == 2);
                        assert(cover_number[msg.outgoing_edge_indices[i][1]] == 2);
                        (*msg_it)[1 + msg.incoming_edge_indices.size() + 2*i] -= omega * 1.0/cover_number[msg.outgoing_edge_indices[i][0]] * reduced_costs[msg.outgoing_edge_indices[i][0]];
                        (*msg_it)[1 + msg.incoming_edge_indices.size() + 2*i + 1] -= omega * 1.0/cover_number[msg.outgoing_edge_indices[i][1]] * reduced_costs[msg.outgoing_edge_indices[i][1]];
                    }
                }

                /*
                vector<REAL> m(1 + msg.incoming_edge_indices.size() + msg.outgoing_edge_indices.size());
                m[0] = 1.0/cover_number[msg.detection_edge] * reduced_costs[msg.detection_edge];
                for(std::size_t i=0; i<msg.incoming_edge_indices.size(); ++i) {
                    m[1+i] = 1.0/cover_number[msg.incoming_edge_indices[i]] * reduced_costs[msg.incoming_edge_indices[i]];
                }
                for(std::size_t i=0; i<msg.outgoing_edge_indices.size(); ++i) {
                    if(msg.outgoing_edge_indices[i][0] == msg.outgoing_edge_indices[i][1]) {
                        m[1 + msg.incoming_edge_indices.size() + i] = 1.0/cover_number[msg.outgoing_edge_indices[i][0]] * reduced_costs[msg.outgoing_edge_indices[i][0]];
                    } else {
                        m[1 + msg.incoming_edge_indices.size() + i] = 1.0/cover_number[msg.outgoing_edge_indices[i][0]] * reduced_costs[msg.outgoing_edge_indices[i][0]] 
                            + 1.0/cover_number[msg.outgoing_edge_indices[i][1]] * reduced_costs[msg.outgoing_edge_indices[i][1]];
                    }
                }

                for(auto& x : m) { x *= omega; }

                (*msg_it) -= m;
                */
            }
        }
    }

    template<typename SOLVER, typename LEFT_FACTOR, typename RIGHT_FACTOR>
    void construct_constraints(SOLVER& s, 
            LEFT_FACTOR& l, typename SOLVER::variable left_detection_var, typename SOLVER::vector left_incoming_vars, typename SOLVER::vector left_outgoing_vars,
            RIGHT_FACTOR& r) const
    {}

private:
    const std::size_t detection_edge;
    vector<std::size_t> incoming_edge_indices;
    vector<std::array<std::size_t,2>> outgoing_edge_indices; // some edges may be division edges, hence there are two of them
};

template<typename cell_detection_constructor_base, typename CELL_DETECTION_FLOW_FACTOR_CONTAINER, typename CELL_DETECTION_FLOW_MESSAGE_CONTAINER>
class cell_detection_flow_constructor : public cell_detection_constructor_base {
public:
    using cell_detection_constructor_base::cell_detection_constructor_base;
    using detection_factor_container = typename cell_detection_constructor_base::detection_factor_container;
    using FMC = typename cell_detection_constructor_base::FMC;

    virtual detection_factor_container* add_detection_hypothesis_impl(LP<FMC>& lp,
            const INDEX timestep, const INDEX hypothesis_id, 
            const REAL detection_cost, const REAL appearance_cost, const REAL disappearance_cost, 
            const INDEX no_incoming_transition_edges, const INDEX no_incoming_division_edges, 
            const INDEX no_outgoing_transition_edges, const INDEX no_outgoing_division_edges)
    {
           const auto cell_detection_counter = this->cumulative_sum_cell_detection_factors[timestep] + hypothesis_id;
           const auto edges = mcf_factor_->GetFactor()->add_detection_hypothesis(cell_detection_counter, 0.0*detection_cost, 0.0*appearance_cost, 0.0*disappearance_cost);

           if(outgoing_edges.size() <= timestep) {
               outgoing_edges.resize(timestep+1);
               incoming_edges.resize(timestep+1);
               dis_appearance_edges.resize(timestep+1);
           } 
           if(outgoing_edges[timestep].size() <= hypothesis_id) {
               outgoing_edges[timestep].resize(hypothesis_id+1);
               incoming_edges[timestep].resize(hypothesis_id+1);
               dis_appearance_edges[timestep].resize(hypothesis_id+1);
           }

           dis_appearance_edges[timestep][hypothesis_id] = edges;

           mcf_test->add_edge(0, 2*cell_detection_counter + 1, 0, 1, appearance_cost);
           mcf_test->add_edge(2*cell_detection_counter + 1, 2*cell_detection_counter + 2, 0, 1, detection_cost);
           mcf_test->add_edge(2*cell_detection_counter + 2, mcf_test->no_nodes()-1, 0, 1, disappearance_cost);

           return cell_detection_constructor_base::add_detection_hypothesis_impl(lp, timestep, hypothesis_id, detection_cost, appearance_cost, disappearance_cost, no_incoming_transition_edges, no_incoming_division_edges, no_outgoing_transition_edges, no_outgoing_division_edges);
    }

    virtual void add_cell_transition_impl(LP<FMC>& lp,
            const INDEX timestep_prev, const INDEX prev_cell, const INDEX timestep_next, const INDEX next_cell, const REAL cost,
            const INDEX outgoing_edge_index, const INDEX no_outgoing_transition_edges, const INDEX no_outgoing_division_edges, 
            const INDEX incoming_edge_index, const INDEX no_incoming_transition_edges, const INDEX no_incoming_division_edges,
            detection_factor_container* out_cell_factor, detection_factor_container* in_cell_factor)
    {
           const auto cell_detection_counter_out = this->cumulative_sum_cell_detection_factors[timestep_prev] + prev_cell;
           const auto cell_detection_counter_in = this->cumulative_sum_cell_detection_factors[timestep_next] + next_cell;
           const std::size_t edge_no = mcf_factor_->GetFactor()->add_cell_transition(cell_detection_counter_out, cell_detection_counter_in, 0.0*cost);
           outgoing_edges[timestep_prev][prev_cell].push_back( {edge_no, edge_no} );
           assert(outgoing_edges[timestep_prev][prev_cell].size() == outgoing_edge_index+1);
           incoming_edges[timestep_next][next_cell].push_back(edge_no);
           assert(incoming_edges[timestep_next][next_cell].size() == incoming_edge_index+1);

           mcf_test->add_edge(2*cell_detection_counter_out + 2, 2*cell_detection_counter_in + 1, 0, 1, cost);

           cell_detection_constructor_base::add_cell_transition_impl(lp,  timestep_prev,  prev_cell,  timestep_next,  next_cell,  cost,  outgoing_edge_index,  no_outgoing_transition_edges,  no_outgoing_division_edges,  incoming_edge_index,  no_incoming_transition_edges,  no_incoming_division_edges, out_cell_factor, in_cell_factor);
    }

    virtual void add_cell_division_impl(LP<FMC>& lp,
            const INDEX timestep_prev, const INDEX prev_cell, const INDEX timestep_next_1, const INDEX next_cell_1, const INDEX timestep_next_2, const INDEX next_cell_2, const REAL cost,
            const INDEX outgoing_edge_index, const INDEX no_outgoing_transition_edges, const INDEX no_outgoing_division_edges, 
            const INDEX incoming_edge_index_1, const INDEX no_incoming_transition_edges_1, const INDEX no_incoming_division_edges_1, 
            const INDEX incoming_edge_index_2, const INDEX no_incoming_transition_edges_2, const INDEX no_incoming_division_edges_2, 
            detection_factor_container* out_cell_factor, detection_factor_container* in_cell_factor_1, detection_factor_container* in_cell_factor_2)
    {
        const auto cell_detection_counter_out = this->cumulative_sum_cell_detection_factors[timestep_prev] + prev_cell;
        const auto cell_detection_counter_in_1 = this->cumulative_sum_cell_detection_factors[timestep_next_1] + next_cell_1;
        const auto cell_detection_counter_in_2 = this->cumulative_sum_cell_detection_factors[timestep_next_2] + next_cell_2;
        auto edges = mcf_factor_->GetFactor()->add_cell_division(cell_detection_counter_out, cell_detection_counter_in_1, cell_detection_counter_in_2, 0.0*cost);
        outgoing_edges[timestep_prev][prev_cell].push_back( {edges[0], edges[1]} );
        assert(outgoing_edges[timestep_prev][prev_cell].size() == no_outgoing_transition_edges + outgoing_edge_index+1); 
        incoming_edges[timestep_next_1][next_cell_1].push_back(edges[0]);
        assert(incoming_edges[timestep_next_1][next_cell_1].size() == no_incoming_transition_edges_1 + incoming_edge_index_1+1);
        incoming_edges[timestep_next_2][next_cell_2].push_back(edges[1]);
        assert(incoming_edges[timestep_next_2][next_cell_2].size() == no_incoming_transition_edges_2 + incoming_edge_index_2+1);

        cell_detection_constructor_base::add_cell_division_impl(lp,  timestep_prev,  prev_cell,  timestep_next_1,  next_cell_1,  timestep_next_2,  next_cell_2,  cost,  outgoing_edge_index,  no_outgoing_transition_edges,  no_outgoing_division_edges,  incoming_edge_index_1,  no_incoming_transition_edges_1,  no_incoming_division_edges_1,  incoming_edge_index_2,  no_incoming_transition_edges_2,  no_incoming_division_edges_2, out_cell_factor, in_cell_factor_1, in_cell_factor_2);
    }

    virtual void begin(LP<FMC>& lp, const std::size_t no_cell_detection_hypotheses, const std::size_t no_cell_transitions, const std::size_t no_cell_divisions)
    {
        cell_detection_constructor_base::begin(lp,  no_cell_detection_hypotheses, no_cell_transitions, no_cell_divisions); 
        // create mcf factor
        mcf_factor_ = lp.template add_factor<CELL_DETECTION_FLOW_FACTOR_CONTAINER>(no_cell_detection_hypotheses, no_cell_transitions, no_cell_divisions); 

        mcf_test = new MCF::SSP<int,REAL>(2 + 2*no_cell_detection_hypotheses + no_cell_transitions, 1 + 2*no_cell_detection_hypotheses + 2*no_cell_transitions);
        mcf_test->add_edge(0, mcf_test->no_nodes()-1, 0, no_cell_detection_hypotheses, 0.0);
        mcf_test->add_node_excess(0, no_cell_detection_hypotheses);
        mcf_test->add_node_excess(mcf_test->no_nodes()-1, -int(no_cell_detection_hypotheses)); 
    }

    virtual void end(LP<FMC>& lp)
    {
        cell_detection_constructor_base::end(lp); 

        // connect cell detection factors with mcf factor
        assert(this->detection_factors_.size() == incoming_edges.size() && this->detection_factors_.size() == outgoing_edges.size());
        for(std::size_t t=0; t<this->detection_factors_.size(); ++t) {
            assert(this->detection_factors_[t].size() == incoming_edges[t].size() && this->detection_factors_[t].size() == outgoing_edges[t].size());
            for(std::size_t c=0; c<this->detection_factors_[t].size(); ++c) {
                const auto node = this->cumulative_sum_cell_detection_factors[t] + c;
                auto incoming = incoming_edges[t][c];
                auto outgoing = outgoing_edges[t][c];
                incoming.push_back( dis_appearance_edges[t][c][0] );
                outgoing.push_back( {dis_appearance_edges[t][c][2], dis_appearance_edges[t][c][2]} );
                auto* m = lp.template add_message<CELL_DETECTION_FLOW_MESSAGE_CONTAINER>(this->detection_factors_[t][c], mcf_factor_, dis_appearance_edges[t][c][1], incoming.begin(), incoming.end(), outgoing.begin(), outgoing.end());
                m->send_message_to_right();
            } 
        }
        //const auto test_cost = mcf_test->solve();
        //mcf_test->print_flow();
        //std::cout << "test cost = " << test_cost << "\n";
        //exit(1);
        //const auto lb = mcf_factor_->LowerBound();
        //std::cout << "lower bound of mcf factor = " << lb << "\n";
        //exit(1);
    }

private:
    CELL_DETECTION_FLOW_FACTOR_CONTAINER* mcf_factor_ = nullptr; 
    MCF::SSP<int,REAL>* mcf_test;

    std::vector<std::vector<std::vector<std::size_t>>> incoming_edges; // for each timestep and for each cell detection, record incoming edge indices from mcf
    std::vector<std::vector<std::vector<std::array<std::size_t,2>>>> outgoing_edges; // for each timestep and for each cell detection, record outgoing edge indices from mcf
    std::vector<std::vector<std::array<std::size_t,3>>> dis_appearance_edges; // appearance edge, detection edge, disappearance edge
};

} // namespace LP_MP

#endif // LP_MP_CELL_DETECTION_FLOW_FACTOR_HXX
