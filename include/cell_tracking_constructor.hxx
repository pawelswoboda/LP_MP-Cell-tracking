#ifndef LP_MP_CELL_TRACKING_CONSTRUCTOR_HXX
#define LP_MP_CELL_TRACKING_CONSTRUCTOR_HXX

#include <vector>
#include "cell_tracking_input.h"

//#include "cell_tracking_rounding.hxx"

// do zrobienia: do not include detection hypotheses that have no incoming and no outgoing edges into the LP

namespace LP_MP {

  // temporary structure which counts how many incoming and outgoing edges are already used by messages for building the model
  struct cell_tracking_input_helper {
      cell_tracking_input_helper(const cell_tracking_input& input)
      {
          transition_edges.reserve(input.no_cells());
          division_edges.reserve(input.no_cells());

          for(const auto & d : input.cell_detections) {
              add_detection( d.no_incoming_transition_edges, d.no_incoming_division_edges,
                      d.no_outgoing_transition_edges, d.no_outgoing_division_edges);
          }
      }

      struct edge_count {
          std::size_t incoming_count = 0;
          std::size_t outgoing_count = 0;

          std::size_t no_incoming = 0;
          std::size_t no_outgoing = 0;
      };

      std::vector<edge_count> transition_edges;
      std::vector<edge_count> division_edges;

      void add_detection(
              const std::size_t no_incoming_transition_edges, const std::size_t no_incoming_division_edges,
              const std::size_t no_outgoing_transition_edges, const std::size_t no_outgoing_division_edges
              )
      {
          transition_edges.push_back({});
          division_edges.push_back({});

          assert(transition_edges.back().incoming_count == 0);
          assert(transition_edges.back().outgoing_count == 0);
          assert(division_edges.back().incoming_count == 0);
          assert(division_edges.back().outgoing_count == 0);

          transition_edges.back().no_incoming = no_incoming_transition_edges;
          transition_edges.back().no_outgoing = no_outgoing_transition_edges;

          division_edges.back().no_incoming = no_incoming_division_edges;
          division_edges.back().no_outgoing = no_outgoing_division_edges;
      }
      std::size_t no_incoming_transition_edges(const std::size_t cell_index) const
      {
          assert(cell_index < transition_edges.size() && transition_edges.size() == division_edges.size());
          return transition_edges[cell_index].no_incoming;
      }
      std::size_t next_incoming_transition_edge(const std::size_t cell_index)
      {
          assert(transition_edges[cell_index].incoming_count < no_incoming_transition_edges(cell_index));
          return transition_edges[cell_index].incoming_count++;
      }

      std::size_t no_incoming_division_edges(const std::size_t cell_index) const
      {
          assert(cell_index < transition_edges.size() && transition_edges.size() == division_edges.size());
          return division_edges[cell_index].no_incoming;
      }
      std::size_t next_incoming_division_edge(const std::size_t cell_index)
      {
          assert(division_edges[cell_index].incoming_count < no_incoming_division_edges(cell_index));
          return division_edges[cell_index].incoming_count++;
      }

      std::size_t no_outgoing_transition_edges(const std::size_t cell_index) const
      {
          assert(cell_index < transition_edges.size() && transition_edges.size() == division_edges.size());
          return transition_edges[cell_index].no_outgoing;
      }
      std::size_t next_outgoing_transition_edge(const std::size_t cell_index)
      {
          assert(transition_edges[cell_index].outgoing_count < no_outgoing_transition_edges(cell_index));
          return transition_edges[cell_index].outgoing_count++;
      }

      std::size_t no_outgoing_division_edges(const std::size_t cell_index) const
      {
          assert(cell_index < transition_edges.size() && transition_edges.size() == division_edges.size());
          return division_edges[cell_index].no_outgoing;
      }
      std::size_t next_outgoing_division_edge(const std::size_t cell_index)
      {
          assert(division_edges[cell_index].outgoing_count < no_outgoing_division_edges(cell_index));
          return division_edges[cell_index].outgoing_count++;
      }
  };

  template<typename DETECTION_FACTOR_CONTAINER>
  class basic_cell_tracking_constructor {
  public:
  using detection_factor_container = DETECTION_FACTOR_CONTAINER;
  using FMC = typename detection_factor_container::FMC;

  using CONSTRUCTOR = basic_cell_tracking_constructor<DETECTION_FACTOR_CONTAINER>;

  template<typename SOLVER>
  basic_cell_tracking_constructor(SOLVER& solver) 
  : lp_(&solver.GetLP()) 
  {}

  virtual detection_factor_container* add_detection_hypothesis_impl(LP<FMC>& lp,
      const INDEX cell_index,
      const REAL detection_cost, const REAL appearance_cost, const REAL disappearance_cost, 
      const INDEX no_incoming_transition_edges, const INDEX no_incoming_division_edges, 
      const INDEX no_outgoing_transition_edges, const INDEX no_outgoing_division_edges) = 0;

  virtual void add_cell_transition_impl(LP<FMC>& lp, const REAL cost,
		  const INDEX cell_outgoing, const INDEX outgoing_edge_index, const INDEX no_outgoing_transition_edges, const INDEX no_outgoing_division_edges, 
          const INDEX cell_incoming, const INDEX incoming_edge_index, const INDEX no_incoming_transition_edges, const INDEX no_incoming_division_edges,
		  detection_factor_container* out_cell_factor, detection_factor_container* in_cell_factor) = 0;

  virtual void add_cell_division_impl(LP<FMC>& lp, const REAL cost,
          const INDEX cell_outgoing, const INDEX outgoing_edge_index, const INDEX no_outgoing_transition_edges, const INDEX no_outgoing_division_edges, 
          const INDEX cell_incoming_1, const INDEX incoming_edge_index_1, const INDEX no_incoming_transition_edges_1, const INDEX no_incoming_division_edges_1, 
          const INDEX cell_incoming_2, const INDEX incoming_edge_index_2, const INDEX no_incoming_transition_edges_2, const INDEX no_incoming_division_edges_2, 
          detection_factor_container* out_cell_factor, detection_factor_container* in_cell_factor_1, detection_factor_container* in_cell_factor_2) = 0;
		
  virtual void add_exclusion_constraint_impl(LP<FMC>& lp, std::vector<DETECTION_FACTOR_CONTAINER*> factors) = 0;

  virtual void begin(LP<FMC>& lp, const std::size_t no_cell_detection_hypotheses, const std::size_t no_transitions, const std::size_t no_divisions)
  {
        if(debug()) {
            std::cout << "no cells = " << no_cell_detection_hypotheses << ", no transitions = " << no_transitions << ", no divisions = " << no_divisions << "\n";
        } 
  }

  virtual void end(LP<FMC>& lp)
  {
    if(debug()) { std::cout << "order factors\n"; }

    //for(INDEX t=0; t<detection_factors_.size(); ++t) {
    //  for(INDEX i=0; i<detection_factors_[t].size()-1; ++i) {
    //      lp.put_in_same_partition(detection_factors_[t][i], detection_factors_[t][i+1]);
    //  }
    //}
  }

  void construct_cell_detection_factors(const cell_tracking_input& input, cell_tracking_input_helper& input_helper)
  {
      assert(std::is_sorted(input.cell_detections.begin(), input.cell_detections.end()));

      cell_detections.reserve(input.cell_detections.size());

      for(std::size_t idx=0; idx<input.cell_detections.size(); ++idx) {
          const auto& d = input.cell_detections[idx];
          auto* f = this->add_detection_hypothesis_impl(*lp_,
                   idx, d.detection_cost, d.appearance_cost, d.disappearance_cost,
                  d.no_incoming_transition_edges, d.no_incoming_division_edges, d.no_outgoing_transition_edges, d.no_outgoing_division_edges);

          cell_detections.push_back({d.timestep, d.cell_number, f});

      }

      for(std::size_t i=0; i<cell_detections.size()-1; ++i) {
          lp_->AddFactorRelation(get_detection_factor(i), get_detection_factor(i+1));
      }

      timestep_to_first_cell_index = input.timestep_to_first_cell_index;
  }

  void add_transitions(const cell_tracking_input& input, cell_tracking_input_helper& input_helper)
  {
      for(const auto& t : input.cell_transitions) {

          auto* out_cell_factor = get_detection_factor( t.outgoing_cell );
          const std::size_t outgoing_edge_index  = input_helper.next_outgoing_transition_edge( t.outgoing_cell );
          const std::size_t no_outgoing_transition_edges = input_helper.no_outgoing_transition_edges( t.outgoing_cell );
          const std::size_t no_outgoing_division_edges = input_helper.no_outgoing_division_edges( t.outgoing_cell );
          out_cell_factor->GetFactor()->set_outgoing_transition_cost(no_outgoing_transition_edges, no_outgoing_division_edges, outgoing_edge_index, 0.5*t.cost);

          auto* in_cell_factor = get_detection_factor( t.incoming_cell );
          const std::size_t incoming_edge_index = input_helper.next_incoming_transition_edge( t.incoming_cell );
          const std::size_t no_incoming_transition_edges = input_helper.no_incoming_transition_edges( t.incoming_cell );
          const std::size_t no_incoming_division_edges = input_helper.no_outgoing_division_edges( t.incoming_cell );
          in_cell_factor->GetFactor()->set_incoming_transition_cost(no_incoming_transition_edges, no_incoming_division_edges, incoming_edge_index, 0.5*t.cost);

          this->add_cell_transition_impl(*lp_, t.cost,
                  t.outgoing_cell, outgoing_edge_index, no_outgoing_transition_edges, no_outgoing_division_edges, 
                  t.incoming_cell, incoming_edge_index, no_incoming_transition_edges, no_incoming_division_edges,
                  out_cell_factor, in_cell_factor);

      }
  }

  void add_divisions(const cell_tracking_input& input, cell_tracking_input_helper& input_helper)
  {
      for(const auto& t : input.cell_divisions) {

          auto* out_cell_factor = get_detection_factor( t.outgoing_cell );
          const std::size_t outgoing_edge_index  = input_helper.next_outgoing_transition_edge( t.outgoing_cell );
          const std::size_t no_outgoing_transition_edges = input_helper.no_outgoing_transition_edges( t.outgoing_cell );
          const std::size_t no_outgoing_division_edges = input_helper.no_outgoing_division_edges( t.outgoing_cell );
          out_cell_factor->GetFactor()->set_outgoing_transition_cost(no_outgoing_transition_edges, no_outgoing_division_edges, outgoing_edge_index, t.cost/3.0);

          auto* in_cell_factor_1 = get_detection_factor( t.incoming_cell_1 );
          const std::size_t incoming_edge_index_1 = input_helper.next_incoming_transition_edge( t.incoming_cell_1 );
          const std::size_t no_incoming_transition_edges_1 = input_helper.no_incoming_transition_edges( t.incoming_cell_1 );
          const std::size_t no_incoming_division_edges_1 = input_helper.no_outgoing_division_edges( t.incoming_cell_1 );
          in_cell_factor_1->GetFactor()->set_incoming_transition_cost(no_incoming_transition_edges_1, no_incoming_division_edges_1, incoming_edge_index_1, t.cost/3.0);

          auto* in_cell_factor_2 = get_detection_factor( t.incoming_cell_2 );
          const std::size_t incoming_edge_index_2 = input_helper.next_incoming_transition_edge( t.incoming_cell_2 );
          const std::size_t no_incoming_transition_edges_2 = input_helper.no_incoming_transition_edges( t.incoming_cell_2 );
          const std::size_t no_incoming_division_edges_2 = input_helper.no_outgoing_division_edges( t.incoming_cell_2 );
          in_cell_factor_2->GetFactor()->set_incoming_transition_cost(no_incoming_transition_edges_2, no_incoming_division_edges_2, incoming_edge_index_2, t.cost/3.0);

          this->add_cell_division_impl(*lp_, t.cost,
                  t.outgoing_cell, outgoing_edge_index, no_outgoing_transition_edges, no_outgoing_division_edges, 
                  t.incoming_cell_1, incoming_edge_index_1, no_incoming_transition_edges_1, no_incoming_division_edges_1,
                  t.incoming_cell_2, incoming_edge_index_2, no_incoming_transition_edges_2, no_incoming_division_edges_2,
                  out_cell_factor, in_cell_factor_1, in_cell_factor_2);

      }

      
  }

  void add_exclusion_constraints(const cell_tracking_input& input, cell_tracking_input_helper& input_helper)
  {
      std::vector<DETECTION_FACTOR_CONTAINER*> conflict_factors;

      for(std::size_t c=0; c<input.no_conflicts(); ++c) {

          conflict_factors.clear();
          auto [conflict_begin, conflict_end] = input.get_conflict(c);
          for(auto it=conflict_begin; it!=conflict_end; ++it) {
              conflict_factors.push_back( get_detection_factor(*it) );
          }
          this->add_exclusion_constraint_impl(*lp_, conflict_factors);

      }
  }

  void construct(const cell_tracking_input& input)
  {
      cell_tracking_input_helper input_helper(input);
      construct_cell_detection_factors(input, input_helper);
      add_transitions(input, input_helper);
      add_divisions(input, input_helper);
      add_exclusion_constraints(input, input_helper); 
  }

public:
  // make protected again
  std::vector<std::size_t> cumulative_sum_cell_detection_factors;

protected:

  struct cell_detection {
      std::size_t timestep;
      std::size_t cell_number;
      DETECTION_FACTOR_CONTAINER* f = nullptr;
  };
  std::vector<cell_detection> cell_detections;
  DETECTION_FACTOR_CONTAINER* get_detection_factor(const std::size_t idx) const
  {
      assert(idx < cell_detections.size());
      return cell_detections[idx].f;
  }

  std::vector<std::size_t> timestep_to_first_cell_index; // given timestep, what is the first cell, assumes that cell_detections is sorted, 
  std::size_t no_timesteps() const { return timestep_to_first_cell_index.size(); }
  std::size_t no_cells(const std::size_t timestep) const {
      assert(timestep < no_timesteps());
      if(timestep < no_timesteps()-1) {
          return timestep_to_first_cell_index[timestep+1] - timestep_to_first_cell_index[timestep];
      } else {
          return no_cells() - timestep_to_first_cell_index[timestep];
      }
  }
  std::size_t cell_index(const std::size_t timestep, const std::size_t cell_number) const 
  {
      assert(timestep < no_timesteps());
      assert(cell_number < no_cells(timestep));
      return timestep_to_first_cell_index[timestep] + cell_number; 
  }

  LP<FMC>* lp_;
};

template<typename BASIC_CELL_TRACKING_CONSTRUCTOR, typename AT_MOST_ONE_CELL_FACTOR_CONTAINER, typename AT_MOST_ONE_CELL_MESSAGE_CONTAINER>
class cell_tracking_exclusion_constructor : public BASIC_CELL_TRACKING_CONSTRUCTOR {
public:
  using detection_factor_container = typename BASIC_CELL_TRACKING_CONSTRUCTOR::detection_factor_container;
  using FMC = typename detection_factor_container::FMC;
  using at_most_one_cell_factor_container = AT_MOST_ONE_CELL_FACTOR_CONTAINER;
  using exclusion_factor = AT_MOST_ONE_CELL_FACTOR_CONTAINER;

  using BASIC_CELL_TRACKING_CONSTRUCTOR::BASIC_CELL_TRACKING_CONSTRUCTOR;

  void add_exclusion_constraint_impl(LP<FMC>& lp, const std::vector<detection_factor_container*> factors)
  {
    auto* e = lp.template add_factor<AT_MOST_ONE_CELL_FACTOR_CONTAINER>(factors.size());
    std::size_t msg_idx = 0;
    for(auto* f : factors) {
        lp.template add_message<AT_MOST_ONE_CELL_MESSAGE_CONTAINER>(f, e, msg_idx++);
        lp.template put_in_same_partition(f,e);

        //lp.ForwardPassFactorRelation(e, f);
        //lp.BackwardPassFactorRelation(e, f);

      //if(f != f_last) {
      //  lp.ForwardPassFactorRelation(f,e);
      //} else {
      //  lp.ForwardPassFactorRelation(e,f);
      //}
      //if(f != f_first) {
      //  lp.BackwardPassFactorRelation(f,e);
      //} else {
      //  lp.BackwardPassFactorRelation(e,f);
      //}
    }
  } 
};



template<typename BASIC_CELL_TRACKING_CONSTRUCTOR>
class cell_tracking_constructor : public BASIC_CELL_TRACKING_CONSTRUCTOR {
public:
  using BASIC_CELL_TRACKING_CONSTRUCTOR::BASIC_CELL_TRACKING_CONSTRUCTOR;

  using detection_factor_container = typename BASIC_CELL_TRACKING_CONSTRUCTOR::detection_factor_container;
  using FMC = typename detection_factor_container::FMC;

  virtual detection_factor_container* add_detection_hypothesis_impl(LP<FMC>& lp,
      const INDEX cell_index, const REAL detection_cost, const REAL appearance_cost, const REAL disappearance_cost, 
      const INDEX no_incoming_transition_edges, const INDEX no_incoming_division_edges, 
      const INDEX no_outgoing_transition_edges, const INDEX no_outgoing_division_edges) 
  {
    auto* f = lp.template add_factor<detection_factor_container>(no_incoming_transition_edges, no_incoming_division_edges, no_outgoing_transition_edges, no_outgoing_division_edges, detection_cost, appearance_cost, disappearance_cost);
    return f;
  }
};

template<typename BASIC_CELL_TRACKING_CONSTRUCTOR, typename TRANSITION_MESSAGE_CONTAINER>
class transition_message_cell_tracking_constructor : public BASIC_CELL_TRACKING_CONSTRUCTOR {
public:
  using BASIC_CELL_TRACKING_CONSTRUCTOR::BASIC_CELL_TRACKING_CONSTRUCTOR;
  using transition_message_container = TRANSITION_MESSAGE_CONTAINER;
  using detection_factor_container = typename BASIC_CELL_TRACKING_CONSTRUCTOR::detection_factor_container;
  using FMC = typename detection_factor_container::FMC;

  //void set_number_of_timesteps(const INDEX t)
  //{
  //  assert(detection_factors_.size() == 0);
  //  detection_factors_.resize(t);
  //}

  virtual void add_cell_transition_impl(LP<FMC>& lp, const REAL cost,
		  const INDEX cell_outgoing, const INDEX outgoing_edge_index, const INDEX no_outgoing_transition_edges, const INDEX no_outgoing_division_edges, 
          const INDEX cell_incoming, const INDEX incoming_edge_index, const INDEX no_incoming_transition_edges, const INDEX no_incoming_division_edges,
		  detection_factor_container* out_cell_factor, detection_factor_container* in_cell_factor)
  {
      lp.template add_message<TRANSITION_MESSAGE_CONTAINER>( out_cell_factor, in_cell_factor, false, outgoing_edge_index, incoming_edge_index);
  }

  virtual void add_cell_division_impl(LP<FMC>& lp, const REAL cost,
          const INDEX cell_outgoing, const INDEX outgoing_edge_index, const INDEX no_outgoing_transition_edges, const INDEX no_outgoing_division_edges, 
          const INDEX cell_incoming_1, const INDEX incoming_edge_index_1, const INDEX no_incoming_transition_edges_1, const INDEX no_incoming_division_edges_1, 
          const INDEX cell_incoming_2, const INDEX incoming_edge_index_2, const INDEX no_incoming_transition_edges_2, const INDEX no_incoming_division_edges_2, 
          detection_factor_container* out_cell_factor, detection_factor_container* in_cell_factor_1, detection_factor_container* in_cell_factor_2)
  {
	  lp.template add_message<TRANSITION_MESSAGE_CONTAINER>(out_cell_factor, in_cell_factor_1, true, no_outgoing_transition_edges + outgoing_edge_index, no_incoming_transition_edges_1 + incoming_edge_index_1);

	  lp.template add_message<TRANSITION_MESSAGE_CONTAINER>(out_cell_factor, in_cell_factor_2, true, no_outgoing_transition_edges + outgoing_edge_index, no_incoming_transition_edges_2 + incoming_edge_index_2);
  }
};

template<typename BASIC_CELL_TRACKING_CONSTRUCTOR, 
  typename MAPPING_EDGE_FACTOR_CONTAINER, typename DIVISION_EDGE_FACTOR_CONTAINER, 
  typename INCOMING_MAPPING_EDGE_MESSAGE_CONTAINER, typename OUTGOING_MAPPING_EDGE_MESSAGE_CONTAINER,
  typename INCOMING_DIVISION_EDGE_MESSAGE_CONTAINER, typename OUTGOING_DIVISION_EDGE_MESSAGE_CONTAINER>
class cell_tracking_constructor_duplicate_edges : public BASIC_CELL_TRACKING_CONSTRUCTOR {
public:
  using CONSTRUCTOR = cell_tracking_constructor_duplicate_edges<BASIC_CELL_TRACKING_CONSTRUCTOR, 
        MAPPING_EDGE_FACTOR_CONTAINER, DIVISION_EDGE_FACTOR_CONTAINER, 
        INCOMING_MAPPING_EDGE_MESSAGE_CONTAINER, OUTGOING_MAPPING_EDGE_MESSAGE_CONTAINER,
        INCOMING_DIVISION_EDGE_MESSAGE_CONTAINER, OUTGOING_DIVISION_EDGE_MESSAGE_CONTAINER>; 

  using BASIC_CELL_TRACKING_CONSTRUCTOR::BASIC_CELL_TRACKING_CONSTRUCTOR;
  using detection_factor_container = typename BASIC_CELL_TRACKING_CONSTRUCTOR::detection_factor_container;
  using FMC = typename detection_factor_container::FMC;
  using mapping_edge_factor_container = MAPPING_EDGE_FACTOR_CONTAINER;
  using division_edge_factor_container = DIVISION_EDGE_FACTOR_CONTAINER;
  using incoming_mapping_edge_message_container = INCOMING_MAPPING_EDGE_MESSAGE_CONTAINER;
  using outgoing_mapping_edge_message_container = OUTGOING_MAPPING_EDGE_MESSAGE_CONTAINER;
  using incoming_division_edge_message_container = INCOMING_DIVISION_EDGE_MESSAGE_CONTAINER;
  using outgoing_division_edge_message_container = OUTGOING_DIVISION_EDGE_MESSAGE_CONTAINER;

/*
  detection_factor_container* add_detection_hypothesis_impl(LP<FMC>& lp,
      const INDEX timestep, const INDEX hypothesis_id,
      const REAL detection_cost, const REAL appearance_cost, const REAL disappearance_cost,
      const INDEX no_incoming_transition_edges, const INDEX no_incoming_division_edges,
      const INDEX no_outgoing_transition_edges, const INDEX no_outgoing_division_edges)
  {
    return new detection_factor_container(no_incoming_transition_edges, no_incoming_division_edges, no_outgoing_transition_edges, no_outgoing_division_edges, detection_cost, appearance_cost, disappearance_cost);
  } 
*/

  void add_cell_transition_impl(LP<FMC>& lp,
		  const INDEX timestep_prev, const INDEX prev_cell, const INDEX timestep_next, const INDEX next_cell, const REAL cost,
		  const INDEX outgoing_edge_index, const INDEX no_incoming_transition_edges, const INDEX no_incoming_division_edges, 
      const INDEX incoming_edge_index, const INDEX no_outgoing_transition_edges, const INDEX no_outgoing_division_edges,
		  detection_factor_container* out_cell_factor, detection_factor_container* in_cell_factor)
  {
    auto* f = lp.template add_factor<MAPPING_EDGE_FACTOR_CONTAINER>(0.0);

    lp.template add_message<OUTGOING_MAPPING_EDGE_MESSAGE_CONTAINER>(f, out_cell_factor, outgoing_edge_index, false);
    lp.template add_message<INCOMING_MAPPING_EDGE_MESSAGE_CONTAINER>(f, in_cell_factor, incoming_edge_index);
    lp.AddFactorRelation(out_cell_factor, f);
    lp.AddFactorRelation(f, in_cell_factor);
  }

  void add_cell_division_impl(LP<FMC>& lp,
		  const INDEX timestep_prev, const INDEX prev_cell, const INDEX timestep_next_1, const INDEX next_cell_1, const INDEX timestep_next_2, const INDEX next_cell_2, const REAL cost,
		  const INDEX outgoing_edge_index, const INDEX no_outgoing_transition_edges, const INDEX no_outgoing_division_edges, 
      const INDEX incoming_edge_index_1, const INDEX no_incoming_transition_edges_1, const INDEX no_incoming_division_edges_1, 
      const INDEX incoming_edge_index_2, const INDEX no_incoming_transition_edges_2, const INDEX no_incoming_division_edges_2, 
		  detection_factor_container* out_cell_factor, detection_factor_container* in_cell_factor_1, detection_factor_container* in_cell_factor_2)
  {
    auto* f_1 = lp.template add_factor<DIVISION_EDGE_FACTOR_CONTAINER>(0.0);
    auto* f_2 = lp.template add_factor<DIVISION_EDGE_FACTOR_CONTAINER>(0.0);

    lp.template add_message<OUTGOING_DIVISION_EDGE_MESSAGE_CONTAINER>(f_1, out_cell_factor, no_outgoing_transition_edges + outgoing_edge_index, true);
    lp.template add_message<OUTGOING_DIVISION_EDGE_MESSAGE_CONTAINER>(f_2, out_cell_factor, no_outgoing_transition_edges + outgoing_edge_index, true);

    lp.template add_message<INCOMING_DIVISION_EDGE_MESSAGE_CONTAINER>(f_1, in_cell_factor_1, no_incoming_transition_edges_1 + incoming_edge_index_1);
    lp.template add_message<INCOMING_DIVISION_EDGE_MESSAGE_CONTAINER>(f_2, in_cell_factor_2, no_incoming_transition_edges_2 + incoming_edge_index_2);

    lp.AddFactorRelation(out_cell_factor, f_1);
    lp.AddFactorRelation(out_cell_factor, f_2);
    lp.AddFactorRelation(f_1, in_cell_factor_1);
    lp.AddFactorRelation(f_2, in_cell_factor_2);
  }
};

template<typename CELL_TRACKING_CONSTRUCTOR>
class cell_tracking_with_division_distance_constructor : public CELL_TRACKING_CONSTRUCTOR {
public:
  using detection_factor_container = typename CELL_TRACKING_CONSTRUCTOR::detection_factor_container;
  using FMC = typename detection_factor_container::FMC;
  //using transition_message_container = typename CELL_TRACKING_CONSTRUCTOR::transition_message_container;
  using CELL_TRACKING_CONSTRUCTOR::CELL_TRACKING_CONSTRUCTOR;

  void set_division_distance(const INDEX d) 
  {
    division_distance_ = d;
  }
  INDEX division_distance() const { return division_distance_; }

  virtual detection_factor_container* add_detection_hypothesis_impl(LP<FMC>& lp,
      const INDEX cell_index, const REAL detection_cost, const REAL appearance_cost, const REAL disappearance_cost,
      const INDEX no_incoming_transition_edges, const INDEX no_incoming_division_edges,
      const INDEX no_outgoing_transition_edges, const INDEX no_outgoing_division_edges)
  {
    assert(division_distance_ >= 2);

    auto* f = lp.template add_factor<detection_factor_container>(no_incoming_transition_edges, no_incoming_division_edges, no_outgoing_transition_edges, no_outgoing_division_edges, detection_cost, appearance_cost, disappearance_cost, division_distance_);
    return f;
  }
  

  /*
  template<typename LP_TYPE>
  void add_cell_division(LP_TYPE& lp, const INDEX timestep_prev, const INDEX prev_cell, const INDEX timestep_next_1, const INDEX next_cell_1, const INDEX timestep_next_2, const INDEX next_cell_2, const REAL cost) 
  {
    auto* out_cell_factor = this->detection_factors_[timestep_prev][prev_cell];
    const INDEX outgoing_edge_index  = this->tc_.next_outgoing_division_edge(timestep_prev, prev_cell);//tc.current_division_no[timestep_prev][prev_cell][1];
    out_cell_factor->GetFactor()->set_outgoing_division_cost(outgoing_edge_index, 1.0/3.0*cost);
    //tc.current_division_no[timestep_prev][prev_cell][1]++;

    auto* in_cell_factor_1 = this->detection_factors_[timestep_next_1][next_cell_1];
    const INDEX incoming_edge_index_1 = this->tc_.next_incoming_division_edge(timestep_next_1, next_cell_1);//tc.current_division_no[timestep_next_1][next_cell_1][0];
    in_cell_factor_1->GetFactor()->set_incoming_division_cost(incoming_edge_index_1, 1.0/3.0*cost);
    //tc.current_division_no[timestep_next_1][next_cell_1][0]++;
    
    auto* in_cell_factor_2 = this->detection_factors_[timestep_next_2][next_cell_2];
    const INDEX incoming_edge_index_2 = this->tc_.next_incoming_division_edge(timestep_next_2, next_cell_2);//tc.current_division_no[timestep_next_2][next_cell_2][0];
    in_cell_factor_2->GetFactor()->set_incoming_division_cost(incoming_edge_index_2, 1.0/3.0*cost);
    //tc.current_division_no[timestep_next_2][next_cell_2][0]++;
    
    auto* m1 = new transition_message_container(out_cell_factor, in_cell_factor_1, true, outgoing_edge_index, incoming_edge_index_1);
    lp.AddMessage(m1);
    
    auto* m2 = new transition_message_container(out_cell_factor, in_cell_factor_2, true, outgoing_edge_index, incoming_edge_index_2);
    lp.AddMessage(m2);

    //std::cout << "DA: " << timestep << " " << prev_cell << ", " << next_cell_1 << " " << next_cell_2 << " " << cost << std::endl;
  }
  */

private:
  INDEX division_distance_ = 0;
};








template<typename CELL_TRACKING_CONSTRUCTOR,
  typename MAPPING_EDGE_FACTOR_CONTAINER, typename DIVISION_EDGE_FACTOR_CONTAINER, 
  typename INCOMING_MAPPING_EDGE_MESSAGE_CONTAINER, typename OUTGOING_MAPPING_EDGE_MESSAGE_CONTAINER,
  typename INCOMING_DIVISION_EDGE_MESSAGE_CONTAINER, typename OUTGOING_DIVISION_EDGE_MESSAGE_CONTAINER>
class cell_tracking_with_division_distance_and_duplicate_edges_constructor : public CELL_TRACKING_CONSTRUCTOR {
public:
  using detection_factor_container = typename CELL_TRACKING_CONSTRUCTOR::detection_factor_container;
  using FMC = typename detection_factor_container::FMC;
  using CELL_TRACKING_CONSTRUCTOR::CELL_TRACKING_CONSTRUCTOR;
  using mapping_edge_factor_container = MAPPING_EDGE_FACTOR_CONTAINER;
  using division_edge_factor_container = DIVISION_EDGE_FACTOR_CONTAINER;
  using incoming_mapping_edge_message_container = INCOMING_MAPPING_EDGE_MESSAGE_CONTAINER;
  using outgoing_mapping_edge_message_container = OUTGOING_MAPPING_EDGE_MESSAGE_CONTAINER;
  using incoming_division_edge_message_container = INCOMING_DIVISION_EDGE_MESSAGE_CONTAINER;
  using outgoing_division_edge_message_container = OUTGOING_DIVISION_EDGE_MESSAGE_CONTAINER;

  void add_cell_transition_impl(LP<FMC>& lp,
                  const INDEX timestep_prev, const INDEX prev_cell, const INDEX timestep_next, const INDEX next_cell, const REAL cost,
                  const INDEX outgoing_edge_index, const INDEX no_outgoing_transition_edges, const INDEX no_outgoing_division_edges, 
                  const INDEX incoming_edge_index, const INDEX no_incoming_transition_edges, const INDEX no_incoming_division_edges,
                  detection_factor_container* out_cell_factor, detection_factor_container* in_cell_factor)
  {
    auto* f = lp.template add_factor<mapping_edge_factor_container>(0.0, this->division_distance());

    lp.template add_message<outgoing_mapping_edge_message_container>(f, out_cell_factor, outgoing_edge_index);
    lp.template add_message<incoming_mapping_edge_message_container>(f, in_cell_factor, incoming_edge_index);

    lp.AddFactorRelation(out_cell_factor, f);
    lp.AddFactorRelation(f, in_cell_factor);
  }

  void add_cell_division_impl(LP<FMC>& lp,
                  const INDEX timestep_prev, const INDEX prev_cell, const INDEX timestep_next_1, const INDEX next_cell_1, const INDEX timestep_next_2, const INDEX next_cell_2, const REAL cost,
                  const INDEX outgoing_edge_index, const INDEX no_outgoing_transition_edges, const INDEX no_outgoing_division_edges, 
                  const INDEX incoming_edge_index_1, const INDEX no_incoming_transition_edges_1, const INDEX no_incoming_division_edges_1, 
                  const INDEX incoming_edge_index_2, const INDEX no_incoming_transition_edges_2, const INDEX no_incoming_division_edges_2, 
                  detection_factor_container* out_cell_factor, detection_factor_container* in_cell_factor_1, detection_factor_container* in_cell_factor_2)
  {
    auto* f_1 = lp.template add_factor<division_edge_factor_container>(0.0);
    auto* f_2 = lp.template add_factor<division_edge_factor_container>(0.0);

    lp.template add_message<outgoing_division_edge_message_container>(f_1, out_cell_factor, outgoing_edge_index);
    lp.template add_message<outgoing_division_edge_message_container>(f_2, out_cell_factor, outgoing_edge_index);

    lp.template add_message<incoming_division_edge_message_container>(f_1, in_cell_factor_1, incoming_edge_index_1);
    lp.template add_message<incoming_division_edge_message_container>(f_2, in_cell_factor_2, incoming_edge_index_2);

    lp.AddFactorRelation(out_cell_factor, f_1);
    lp.AddFactorRelation(out_cell_factor, f_2);
    lp.AddFactorRelation(f_1, in_cell_factor_1);
    lp.AddFactorRelation(f_2, in_cell_factor_2);
  } 
};








// this constructor first builds an ordinary cell tracking problem, solves it, and afterwards converts it into a cell tracking with division distance problem and solves it when it is preprocessed
// version with no duplicate edges
template<
typename CELL_TRACKING_CONSTRUCTOR,
typename CELL_TRACKING_DIVISION_DISTANCE_CONSTRUCTOR >
class cell_tracking_division_distance_conversion_constructor
{
public:
  using detection_factor_container = typename CELL_TRACKING_CONSTRUCTOR::detection_factor_container;
  using FMC = typename detection_factor_container::FMC;
  using at_most_one_cell_factor_container = typename CELL_TRACKING_CONSTRUCTOR::at_most_one_cell_factor_container;
  using exclusion_factor = typename CELL_TRACKING_CONSTRUCTOR::exclusion_factor;

  cell_tracking_division_distance_conversion_constructor()
    : cdc_()
  {}

  ~cell_tracking_division_distance_conversion_constructor()
  {
    if(ctc_dd_ != nullptr) {
      delete ctc_dd_;
    }
  }

  template<typename LP_TYPE>
  detection_factor_container* 
  add_detection_hypothesis(
      LP_TYPE& lp, 
      const INDEX timestep, const INDEX hypothesis_id, 
      const REAL detection_cost, const REAL appearance_cost, const REAL disappearance_cost, 
      const INDEX no_incoming_transition_edges, const INDEX no_incoming_division_edges, 
      const INDEX no_outgoing_transition_edges, const INDEX no_outgoing_division_edges
      )
  { 
    if(timestep >= no_transition_edges_.size()) {
      no_transition_edges_.resize(timestep+1);
    }
    if(hypothesis_id >= no_transition_edges_[timestep].size()) {
      no_transition_edges_[timestep].resize(hypothesis_id + 1, {0,0});
    }
    no_transition_edges_[timestep][hypothesis_id][0] = no_incoming_transition_edges;
    no_transition_edges_[timestep][hypothesis_id][1] = no_outgoing_transition_edges;

    cdc_.add_detection_hypothesis(
        lp, timestep, hypothesis_id, detection_cost, appearance_cost, disappearance_cost, 
        no_incoming_transition_edges, no_incoming_division_edges,
        no_outgoing_transition_edges, no_outgoing_division_edges
        );
  }

  template<typename LP_TYPE>
  void add_cell_transition(LP_TYPE& lp, const INDEX timestep_prev, const INDEX prev_cell, const INDEX timestep_next, const INDEX next_cell, const REAL cost) 
  {
    cdc_.add_cell_transition.push_back(lp, timestep_prev, prev_cell, timestep_next, next_cell, cost);
    transition_edges_.push_back({timestep_prev, prev_cell, timestep_next, next_cell});

  }

  template<typename LP_TYPE>
  void add_cell_division(LP_TYPE& lp, const INDEX timestep_prev, const INDEX prev_cell, const INDEX timestep_next_1, const INDEX next_cell_1, const INDEX timestep_next_2, const INDEX next_cell_2, const REAL cost) 
  {
    cdc_.add_cell_division(lp, timestep_prev, prev_cell, timestep_next_1, next_cell_1, timestep_next_2, next_cell_2, cost);
    assert(timestep_prev + 1 == timestep_next_1);
    assert(timestep_prev + 1 == timestep_next_2);
    division_edges_.push_back({timestep_prev, prev_cell, timestep_next_1, next_cell_1, timestep_next_2, next_cell_2});
  }

  template<typename LP_TYPE, typename ITERATOR>
  at_most_one_cell_factor_container* add_exclusion_constraint(LP_TYPE& lp, ITERATOR begin, ITERATOR end) // iterator points to std::array<INDEX,2>
  {
    auto* f = cdc_.add_exclusion_constraint(lp, begin, end);
    exclusion_factors_.push_back(f);
  }

  template<typename LP>
  void convert(LP* lp)
  {
    auto* ctc_dd = new CELL_TRACKING_DIVISION_DISTANCE_CONSTRUCTOR();
    // read all factors from ordinary cell tracking problem, and add equivalent ones to problem with division distance.
    // cost of new factors is reparametrized cost of old ones

    // detection factors
    /*
    for(INDEX t=0; t<cdc_.detection_factors_.size(); ++t) {
      for(INDEX i=0; i<cdc_.detection_factors_[i].size(); ++i) {
        auto* fp = cdc_.detection_factors_[i][t];
        auto& f = *(fp->GetFactor());
        const REAL detection_cost = f.detection_cost();
        const REAL appearance_cost = f.appearance_cost();
        const REAL disappearance_cost = f.disappearance_cost();
        const INDEX no_incoming_transition = no_transition_edges_[i][t][0];
        const INDEX no_outgoing_transition = no_transition_edges_[i][t][1];
        const INDEX no_incoming_division = f.no_incoming_edges() - 1 - no_incoming_transition;
        const INDEX no_outgoing_division = f.no_outgoing_edges() - 1 - no_outgoing_transition;

        auto* fp_dd = ctc_dd_->add_detection_hypothesis( 
            lp, t, 
            detection_cost, appearance_cost, disappearance_cost,
            no_incoming_transition, no_incoming_division,
            no_outgoing_transition, no_outgoing_division
            );
        auto& f_dd = *(fp_dd->GetFactor());

        for(INDEX incoming_edge = 0; incoming_edge<no_incoming_transition; ++incoming_edge) {
          const REAL cost = f.incoming[incoming_edge];
          f_dd.set_incoming_transition_cost(incoming_edge, cost); 
        }
        for(INDEX incoming_edge = 0; incoming_edge<no_incoming_division; ++incoming_edge) {
          const REAL cost = f.incoming[ incoming_edge + no_transition_edges_[t][i] ] ;
          f_dd.set_incoming_division_cost(incoming_edge, cost); 
        }

        for(INDEX outgoing_edge = 0; outgoing_edge<no_outgoing_transition; ++outgoing_edge) {
          const REAL cost = f.outgoing[outgoing_edge];
          f_dd.set_outgoing_transition_cost(outgoing_edge, cost); 
        }
        for(INDEX outgoing_edge = 0; outgoing_edge<no_outgoing_division; ++outgoing_edge) {
          const REAL cost = f.outgoing[ outgoing_edge + no_transition_edges_[t][i] ] ;
          f_dd.set_outgoing_division_cost(outgoing_edge, cost); 
        } 
      } 
    }

    // link detections via edges
    {
      std::vector<std::vector<std::array<INDEX,2>>> edge_counter(cdc_.detection_factors_.size());
      for(INDEX i=0; i<edge_counter.size(); ++i) {
        edge_counter[i].resize( cdc_.detection_factors_[i].size(), {0,0} );
      }
      for(const auto t : transition_edges_) {
        auto* f_prev = ctc_dd->detection_factors_[t[0]][t[1]];
        auto* f_next = ctc_dd_->detection_factors_[t[2]][t[3]];
        const INDEX outgoing_edge_index = edge_counter[t[0]][t[1]][1]++;
        const INDEX incoming_edge_index = edge_counter[t[2]][t[3]][0]++;

        auto* m = new transition_message_container(f_prev, f_next, false, outgoing_edge_index, incoming_edge_index);
        lp.AddMessage(m);
      }
      for(const auto t : division_edges) {
        auto* f_prev = ctc_dd_->detection_factors_[t[0]][t[1]];
        auto* f_next_1 = ctc_dd->detection_factors_[t[2]][t[3]];
        auto* f_next_2 = ctc_dd->detection_factors_[t[4]][t[5]];
        const INDEX outgoing_edge_index = edge_counter[t[0]][t[1]][1]++;
        const INDEX incoming_edge_index_1 = edge_counter[t[2]][t[3]][0]++;
        const INDEX incoming_edge_index_2 = edge_counter[t[4]][t[5]][0]++;

        auto* m1 = new transition_message_container(f_prev, f_next_1, true, outgoing_edge_index, incoming_edge_index_1);
        lp.AddMessage(m1);

        auto* m2 = new transition_message_container(f_prev, f_next_2, true, outgoing_edge_index, incoming_edge_index_2);
        lp.AddMessage(m2); 
      }
    }
    */

    // exclusion factors
    /*
    {
      auto ef_it = exclusion_factors_.begin();
      for(INDEX i=0; i<cdc_.exclusions_.size();) {
        auto item_begin = cdc_.exclusions_.begin() + i;
        auto item_end = item_begin;
        while((*item_end) != exclusion_item_delimiter) {
          ++item_end; 
        }
        auto* f = ctc_dd_->add_exclusion_constraint(*lp_, item_begin, item_end); 
        i += 1 + std::distance(item_begin, item_end);

        auto* f_old = *ef_it;
        ++ef_it;
        for(INDEX j=0; j<f_old.size(); ++j) {
          f[j] = f_old[j];
        } 
      }
      assert(ef_it == exclusions_factors_.end());
    }
    */
    
    // free up space taken up by conversion information
    std::swap(no_transition_edges_,decltype(no_transition_edges_){});
    std::swap(transition_edges_,decltype(transition_edges_){});
    std::swap(division_edges_,decltype(division_edges_){});
    std::swap(exclusion_factors_,decltype(exclusion_factors_){});

  }
protected:
  // store the number of incoming/outgoing transition edges here
  std::vector<std::vector<std::array<INDEX,2>>> no_transition_edges_; 
  std::vector<std::array<INDEX,4>> transition_edges_;
  std::vector<std::array<INDEX,6>> division_edges_;

  std::vector<exclusion_factor*> exclusion_factors_;

  CELL_TRACKING_CONSTRUCTOR cdc_;
  CELL_TRACKING_DIVISION_DISTANCE_CONSTRUCTOR* ctc_dd_ = nullptr;
};


template<typename BASIC_CELL_TRACKING_CONSTRUCTOR, typename DETECTION_MESSAGE_CONTAINER, typename TRANSITION_MESSAGE_CONTAINER, typename AT_MOST_ONE_CELL_FACTOR_CONTAINER, typename AT_MOST_ONE_CELL_MESSAGE_CONTAINER_INCOMING, typename AT_MOST_ONE_CELL_MESSAGE_CONTAINER_OUTGOING>
class cell_tracking_fine_decomposition_constructor : public BASIC_CELL_TRACKING_CONSTRUCTOR {
public:
  using BASIC_CELL_TRACKING_CONSTRUCTOR::BASIC_CELL_TRACKING_CONSTRUCTOR;
  using detection_factor_container = typename BASIC_CELL_TRACKING_CONSTRUCTOR::detection_factor_container;
  using FMC = typename detection_factor_container::FMC;

  virtual detection_factor_container* add_detection_hypothesis_impl(LP<FMC>& lp,
      const INDEX cell_index,
      const REAL detection_cost, const REAL appearance_cost, const REAL disappearance_cost, 
      const INDEX no_incoming_transition_edges, const INDEX no_incoming_division_edges, 
      const INDEX no_outgoing_transition_edges, const INDEX no_outgoing_division_edges) 
  {
    using incoming_container_type = typename detection_factor_container::incoming_factor_type;
    const INDEX incoming_size = no_incoming_transition_edges + no_incoming_division_edges + 1;
    auto* incoming_container = lp.template add_factor<incoming_container_type>(incoming_size);
    auto& incoming_factor = *incoming_container->GetFactor();
    incoming_factor.edges[incoming_size-1] = appearance_cost;
    incoming_factor.detection = 0.5*detection_cost;

    using outgoing_container_type = typename detection_factor_container::outgoing_factor_type;
    const INDEX outgoing_size = no_outgoing_transition_edges + no_outgoing_division_edges + 1;
    auto* outgoing_container = lp.template add_factor<outgoing_container_type>(outgoing_size); 
    auto& outgoing_factor = *outgoing_container->GetFactor();
    outgoing_factor.edges[outgoing_size-1] = disappearance_cost;
    outgoing_factor.detection = 0.5*detection_cost;

    lp.template add_message<DETECTION_MESSAGE_CONTAINER>(incoming_container, outgoing_container);

    return new detection_factor_container({incoming_container, outgoing_container});
  }

  virtual void add_exclusion_constraint_impl(LP<FMC>& lp, const std::vector<detection_factor_container*> factors) // iterator points to std::array<INDEX,2>
  {
    auto* e = lp.template add_factor<AT_MOST_ONE_CELL_FACTOR_CONTAINER>(factors.size());
    INDEX msg_idx = 0;
    for(auto* f : factors) {
      auto* f_incoming = f->incoming;
      auto* f_outgoing = f->outgoing;
      
      lp.template add_message<AT_MOST_ONE_CELL_MESSAGE_CONTAINER_INCOMING>(f_incoming, e, msg_idx);
      lp.template add_message<AT_MOST_ONE_CELL_MESSAGE_CONTAINER_OUTGOING>(f_outgoing, e, msg_idx++);
      lp.put_in_same_partition(f,e);

      lp.AddFactorRelation(f_incoming, f_outgoing);

      //if(f != f_last) {
      //  lp.ForwardPassFactorRelation(f,e);
      //} else {
      //  lp.ForwardPassFactorRelation(e,f);
      //}
      //if(f != f_first) {
      //  lp.BackwardPassFactorRelation(f,e);
      //} else {
      //  lp.BackwardPassFactorRelation(e,f);
      //}
    }
  } 

  virtual void add_cell_transition_impl(LP<FMC>& lp,
		  const INDEX timestep_prev, const INDEX prev_cell, const INDEX timestep_next, const INDEX next_cell, const REAL cost,
		  const INDEX outgoing_edge_index, const INDEX no_outgoing_transition_edges, const INDEX no_outgoing_division_edges, 
      const INDEX incoming_edge_index, const INDEX no_incoming_transition_edges, const INDEX no_incoming_division_edges,
		  detection_factor_container* out_cell_factor, detection_factor_container* in_cell_factor)
  {
    auto* f_outgoing = out_cell_factor->outgoing;
    auto* f_incoming = in_cell_factor->incoming;
    lp.template add_message<TRANSITION_MESSAGE_CONTAINER>(f_outgoing, f_incoming, outgoing_edge_index, incoming_edge_index, false);
  }

  virtual void add_cell_division_impl(LP<FMC>& lp,
		  const INDEX timestep_prev, const INDEX prev_cell, const INDEX timestep_next_1, const INDEX next_cell_1, const INDEX timestep_next_2, const INDEX next_cell_2, const REAL cost,
		  const INDEX outgoing_edge_index, const INDEX no_outgoing_transition_edges, const INDEX no_outgoing_division_edges, 
      const INDEX incoming_edge_index_1, const INDEX no_incoming_transition_edges_1, const INDEX no_incoming_division_edges_1, 
      const INDEX incoming_edge_index_2, const INDEX no_incoming_transition_edges_2, const INDEX no_incoming_division_edges_2, 
		  detection_factor_container* out_cell_factor, detection_factor_container* in_cell_factor_1, detection_factor_container* in_cell_factor_2)
  {
    auto* f_outgoing = out_cell_factor->outgoing;
    auto* f_incoming_1 = in_cell_factor_1->incoming;
    auto* f_incoming_2 = in_cell_factor_2->incoming;

    lp.template add_message<TRANSITION_MESSAGE_CONTAINER>(f_outgoing, f_incoming_1, outgoing_edge_index, incoming_edge_index_1, true);
    lp.template add_message<TRANSITION_MESSAGE_CONTAINER>(f_outgoing, f_incoming_2, outgoing_edge_index, incoming_edge_index_2, true);
  } 
};



// we assume that cells are ordered by heigth and lower cells cannot exit if higher ones do not do so.
template<typename CELL_TRACKING_CONSTRUCTOR,
  typename EXIT_CONSTRAINT_FACTOR,
  typename EXIT_CONSTRAINT_LOWER_MESSAGE,
  typename EXIT_CONSTRAINT_UPPER_MESSAGE
  >
class cell_tracking_mother_machine_constructor : public CELL_TRACKING_CONSTRUCTOR {
public:
  using CELL_TRACKING_CONSTRUCTOR::CELL_TRACKING_CONSTRUCTOR;
  using FMC = typename EXIT_CONSTRAINT_FACTOR::FMC;
  template<typename LP_TYPE>
  void add_exit_constraint(LP_TYPE& lp, const INDEX timestep, const INDEX lower_cell_detection, const INDEX upper_cell_detection)
  {
    if(this->detection_factors_[timestep][lower_cell_detection] == nullptr || this->detection_factors_[timestep][upper_cell_detection] == nullptr) { 
      return;
    }
    //std::cout << "t = " << timestep << ", lower cell = " << lower_cell_detection << ", upper cell = " << upper_cell_detection;
    if(this->detection_factors_[timestep][upper_cell_detection]->GetFactor()->outgoing.size() > 1) { // only if there are outgoing edges which are not exit ones do we need exit constraints
      //std::cout << " add factor";
      auto* e = lp.template add_factor<EXIT_CONSTRAINT_FACTOR>();
      lp.template add_message<EXIT_CONSTRAINT_LOWER_MESSAGE>(this->detection_factors_[timestep][lower_cell_detection], e);
      lp.template add_message<EXIT_CONSTRAINT_UPPER_MESSAGE>(this->detection_factors_[timestep][upper_cell_detection], e);
    }
    //std::cout << "\n";
  } 
};










} // end namespace LP_MP

#endif // LP_MP_CELL_TRACKING_CONSTRUCTOR_HXX

