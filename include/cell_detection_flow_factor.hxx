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
        const std::size_t no_mcf_nodes = 2 + no_cells;    
        // additional edges due to (dis)appearances, cell divisions, source to sink
        const std::size_t no_mcf_edges = no_cell_transitions + 2*no_cells + 2*no_cell_divisions + 1;
        mcf_ = MCF::SSP<int,REAL>(no_mcf_nodes, no_mcf_edges); 

        mcf_.add_node_excess(0, no_cells);
        mcf_.add_node_excess(no_mcf_nodes-1, -no_cells);
        mcf_.add_edge(0, no_mcf_nodes-1, 0, no_cells, 0.0);
    }

    std::array<std::size_t,2> add_detection_hypothesis(const std::size_t cell_detection_counter, const double detection_cost, const double appearance_cost, const double disappearance_cost)
    {
        assert(cell_detection_counter < mcf_.no_nodes()-2);
        // appearance edge
        const auto appearance_edge = mcf_.add_edge(0, cell_detection_counter+1 , 0, 1, appearance_cost);
        // disappearance
        const auto disappearance_edge = mcf_.add_edge(cell_detection_counter+1 , mcf_.no_nodes()-1, 0, 1, disappearance_cost); 
        return {appearance_edge, disappearance_edge};
    }

    std::size_t add_cell_transition(const std::size_t cell_detection_counter_out, const std::size_t cell_detection_counter_in, const double cost)
    {
        return mcf_.add_edge(cell_detection_counter_out, cell_detection_counter_in, 0, 1, cost);
    }

    std::array<std::size_t,2> add_cell_division(const std::size_t cell_detection_counter_out, const std::size_t cell_detection_counter_in_1, const std::size_t cell_detection_counter_in_2, const double cost)
    {
        std::array<std::size_t,2> edges;
        edges[0] = mcf_.add_edge(cell_detection_counter_out, cell_detection_counter_in_1, 0, 1, 0.5*cost);
        edges[1] = mcf_.add_edge(0, cell_detection_counter_in_1, 0, 1, 0.5*cost); // copied edge
        return edges;
    } 

    REAL LowerBound() const
    {
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
    cell_detection_flow_message(INCOMING_EDGE_ITERATOR incoming_begin, INCOMING_EDGE_ITERATOR incoming_end, OUTGOING_EDGES_ITERATOR outgoing_begin, OUTGOING_EDGES_ITERATOR outgoing_end)
    : incoming_edge_indices(incoming_begin, incoming_end),
    outgoing_edge_indices(outgoing_begin, outgoing_end)
    {}

    template<typename LEFT_FACTOR, typename MSG>
    void RepamLeft(LEFT_FACTOR& l, const MSG& msg)
    {
        assert(msg.size() == l.incoming.size() + l.outgoing.size());
        for(std::size_t i=0; i<incoming_edge_indices.size(); ++i) {
            l.update_incoming(i, msg[i]);
        } 

        for(std::size_t i=0; i<outgoing_edge_indices.size(); ++i) {
            l.update_outgoing(i, msg[incoming_edge_indices.size() + i]);
        } 
    }

    template<typename RIGHT_FACTOR, typename MSG>
    void RepamRight(RIGHT_FACTOR& r, const MSG& msg)
    {
        // first messages are incoming
        for(std::size_t i=0; i<incoming_edge_indices.size(); ++i) {
            r.mcf_.update_cost(incoming_edge_indices[i], msg[i]);
        } 

        // then come outgoing edges. Note: duplicated outgoing edge has its edge number recorded in outgoing edge indices
        for(std::size_t i=0; i<outgoing_edge_indices.size(); ++i) {
            if(outgoing_edge_indices[i][0] == outgoing_edge_indices[i][1]) {
                r.mcf_.update_cost(outgoing_edge_indices[i][0], msg[i]);
            } else {
                r.mcf_.update_cost(outgoing_edge_indices[i][0], 0.5*msg[i]);
                r.mcf_.update_cost(outgoing_edge_indices[i][1], 0.5*msg[i]);
            }
        }
    }

    template<typename LEFT_FACTOR, typename MSG>
    void send_message_to_right(LEFT_FACTOR& l, MSG& msg, const REAL omega)
    {
        const auto min_incoming_cost = l.min_incoming();
        const auto min_outgoing_cost = l.min_outgoing();
        vector<REAL> m(l.incoming.size() + l.outgoing.size());
        const auto min_cost = l.min_detection_cost();
        if(min_cost >= 0) {
            for(std::size_t i=0; i<l.incoming.size(); ++i) { m[i] = l.incoming[i] + 0.5*l.detection; }
            for(std::size_t i=0; i<l.outgoing.size(); ++i) { m[i + l.incoming.size()] = 0.5*l.detection + l.outgoing[i]; }
        } else {
            // compute uniform minorant
            assert(false);
            for(std::size_t i=0; i<l.incoming.size(); ++i) { m[i] = l.incoming[i] + 0.5*min_cost; }
            for(std::size_t i=0; i<l.outgoing.size(); ++i) { m[i + l.incoming.size()] = l.outgoing[i] + 0.5*min_cost; } 
        }

        for(auto& x : m) { x *= omega; }

        msg -= m;
    }

    template<typename RIGHT_FACTOR, typename MSG>
    void send_message_to_left(RIGHT_FACTOR& r, MSG& msg, const REAL omega)
    {
        assert(false);
    }

    template<typename RIGHT_FACTOR, typename MSG_ITERATOR>
    static void SendMessagesToLeft(RIGHT_FACTOR& r, MSG_ITERATOR msg_begin, MSG_ITERATOR msg_end, const REAL omega)
    {
        r.mcf_.solve();
        for(auto msg_it=msg_begin; msg_it!=msg_end; ++msg_it) {
            auto& msg = (*msg_it).GetMessageOp();
            vector<REAL> m(msg.incoming_edge_indices.size() + msg.outgoing_edge_indices.size());
            for(std::size_t i=0; i<msg.incoming_edge_indices.size(); ++i) {
                m[i] = r.mcf_.reduced_cost(msg.incoming_edge_indices[i]);
            }
            for(std::size_t i=0; i<msg.outgoing_edge_indices.size(); ++i) {
                if(msg.outgoing_edge_indices[i][0] == msg.outgoing_edge_indices[i][1]) {
                    m[msg.incoming_edge_indices.size() + i] = r.mcf_.reduced_cost(msg.outgoing_edge_indices[i][0]);
                } else {
                    m[msg.incoming_edge_indices.size() + i] = r.mcf_.reduced_cost(msg.outgoing_edge_indices[i][0]) + r.mcf_.reduced_cost(msg.outgoing_edge_indices[i][1]);
                }
            }

            for(auto& x : m) { x *= omega; }

            (*msg_it) -= m;
        }
    }

    template<typename SOLVER, typename LEFT_FACTOR, typename RIGHT_FACTOR>
    void construct_constraints(SOLVER& s, 
            LEFT_FACTOR& l, typename SOLVER::variable left_detection_var, typename SOLVER::vector left_incoming_vars, typename SOLVER::vector left_outgoing_vars,
            RIGHT_FACTOR& r) const
    {}

private:
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
      const auto edges = mcf_factor_->GetFactor()->add_detection_hypothesis(cell_detection_counter, detection_cost, 0.5*appearance_cost, 0.5*disappearance_cost);

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

      return cell_detection_constructor_base::add_detection_hypothesis_impl(lp, timestep, hypothesis_id, detection_cost, 0.5*appearance_cost, 0.5*disappearance_cost, no_incoming_transition_edges, no_incoming_division_edges, no_outgoing_transition_edges, no_outgoing_division_edges);
  }

  virtual void add_cell_transition_impl(LP<FMC>& lp,
		  const INDEX timestep_prev, const INDEX prev_cell, const INDEX timestep_next, const INDEX next_cell, const REAL cost,
		  const INDEX outgoing_edge_index, const INDEX no_outgoing_transition_edges, const INDEX no_outgoing_division_edges, 
          const INDEX incoming_edge_index, const INDEX no_incoming_transition_edges, const INDEX no_incoming_division_edges,
		  detection_factor_container* out_cell_factor, detection_factor_container* in_cell_factor)
  {
      const auto cell_detection_counter_out = this->cumulative_sum_cell_detection_factors[timestep_prev] + prev_cell;
      const auto cell_detection_counter_in = this->cumulative_sum_cell_detection_factors[timestep_next] + next_cell;
      const std::size_t edge_no = mcf_factor_->GetFactor()->add_cell_transition(cell_detection_counter_out, cell_detection_counter_in, 0.5*cost);
      outgoing_edges[timestep_prev][prev_cell].push_back( {edge_no, edge_no} );
      incoming_edges[timestep_next][next_cell].push_back(edge_no);

      cell_detection_constructor_base::add_cell_transition_impl(lp,  timestep_prev,  prev_cell,  timestep_next,  next_cell,  0.5*cost,  outgoing_edge_index,  no_outgoing_transition_edges,  no_outgoing_division_edges,  incoming_edge_index,  no_incoming_transition_edges,  no_incoming_division_edges, out_cell_factor, in_cell_factor);
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
      auto edges = mcf_factor_->GetFactor()->add_cell_division(cell_detection_counter_out, cell_detection_counter_in_1, cell_detection_counter_in_2, 0.5*cost);
      outgoing_edges[timestep_prev][prev_cell].push_back( {edges[0], edges[1]} );
      incoming_edges[timestep_next_1][next_cell_1].push_back(edges[0]);
      incoming_edges[timestep_next_2][next_cell_2].push_back(edges[1]);

      cell_detection_constructor_base::add_cell_division_impl(lp,  timestep_prev,  prev_cell,  timestep_next_1,  next_cell_1,  timestep_next_2,  next_cell_2,  cost,  outgoing_edge_index,  no_outgoing_transition_edges,  no_outgoing_division_edges,  incoming_edge_index_1,  no_incoming_transition_edges_1,  no_incoming_division_edges_1,  incoming_edge_index_2,  no_incoming_transition_edges_2,  no_incoming_division_edges_2, out_cell_factor, in_cell_factor_1, in_cell_factor_2);
  }

  virtual void begin(LP<FMC>& lp, const std::size_t no_cell_detection_hypotheses, const std::size_t no_cell_transitions, const std::size_t no_cell_divisions)
  {
      std::cout << "no cells = " << no_cell_detection_hypotheses << ", no transitions = " << no_cell_transitions << ", no divisions = " << no_cell_divisions << "\n";
      cell_detection_constructor_base::begin(lp, no_cell_detection_hypotheses,no_cell_transitions,no_cell_divisions); 
      // create mcf factor
      mcf_factor_ = lp.template add_factor<CELL_DETECTION_FLOW_FACTOR_CONTAINER>(no_cell_detection_hypotheses,no_cell_transitions,no_cell_divisions); 
  }

  virtual void end(LP<FMC>& lp)
  {
      cell_detection_constructor_base::end(lp); 

      // connect cell detection factors with mcf factor
      assert(this->detection_factors_.size() == incoming_edges.size() && this->detection_factors_.size() == outgoing_edges.size());
      for(std::size_t t=0; t<this->detection_factors_.size(); ++t) {
          assert(this->detection_factors_[t].size() == incoming_edges[t].size() && this->detection_factors_[t].size() == outgoing_edges[t].size());
          for(std::size_t c=0; c<this->detection_factors_[c].size(); ++c) {
              auto& incoming = incoming_edges[t][c];
              auto& outgoing = outgoing_edges[t][c];
              incoming.push_back( dis_appearance_edges[t][c][0] );
              outgoing.push_back( {dis_appearance_edges[t][c][1], dis_appearance_edges[t][c][1]} );
              lp.template add_message<CELL_DETECTION_FLOW_MESSAGE_CONTAINER>(this->detection_factors_[t][c], mcf_factor_, incoming.begin(), incoming.end(), outgoing.begin(), outgoing.end());
          } 
      }
  }

private:
  CELL_DETECTION_FLOW_FACTOR_CONTAINER* mcf_factor_ = nullptr; 
  std::vector<std::vector<std::vector<std::size_t>>> incoming_edges; // for each timestep and for each cell detection, record incoming edge indices from mcf
  std::vector<std::vector<std::vector<std::array<std::size_t,2>>>> outgoing_edges; // for each timestep and for each cell detection, record outgoing edge indices from mcf
  std::vector<std::vector<std::array<std::size_t,2>>> dis_appearance_edges;
};

} // namespace LP_MP

#endif // LP_MP_CELL_DETECTION_FLOW_FACTOR_HXX
