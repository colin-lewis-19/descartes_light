/*
 * Software License Agreement (Apache License)
 *
 * Copyright (c) 2021, Southwest Research Institute
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
#ifndef DESCARTES_LIGHT_SOLVERS_BGL_IMPL_ADD_ALL_DYNAMIC_DIJKSTRA_SOLVER_HPP
#define DESCARTES_LIGHT_SOLVERS_BGL_IMPL_ADD_ALL_DYNAMIC_DIJKSTRA_SOLVER_HPP

#include <descartes_light/descartes_macros.h>
DESCARTES_IGNORE_WARNINGS_PUSH
#include <console_bridge/console.h>
#include <omp.h>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <descartes_light/solvers/bgl/dfs_add_all_solver.h>
#include <descartes_light/types.h>

namespace descartes_light
{

template <typename FloatType>
std::vector<VertexDesc<FloatType>> DFSAddAllSolver<FloatType>::reconstructPath(const VertexDesc<FloatType>& source,
                                                                                             const VertexDesc<FloatType>& target) const
{
  // Reconstruct the path from predecessors
  std::vector<VertexDesc<FloatType>> path;

  VertexDesc<FloatType> v = target;
  path.push_back(v);

  for (VertexDesc<FloatType> u = predecessor_map_.at(v); u != v; v = u, u = predecessor_map_.at(v))
  {
    path.push_back(u);
  }
  std::reverse(path.begin(), path.end());

  // Check that the last traversed vertex is the source vertex
  if (v != source)
    throw std::runtime_error("Failed to find path through the graph");

  return path;
}


template <typename FloatType>
class AddAllVisitor : public boost::default_dijkstra_visitor
{
  public:
  AddAllVisitor(std::vector<typename EdgeEvaluator<FloatType>::ConstPtr>& edge_eval,
                std::vector<std::vector<VertexDesc<FloatType>>>& ladder_rungs,
                BGLGraph<FloatType>& eg)
                : eval_(edge_eval),
                  ladder_rungs_(ladder_rungs),
                  mutable_graph_(eg)
       {
       }

  void examine_vertex(VertexDesc<FloatType> u, BGLGraph<FloatType> g)
  {
    int out_deg = boost::out_degree(u, g);
    // return if the vertex has any out edges
    if (out_deg == 0)
    {
    std::size_t current_rung = g[u].rung_idx;
    std::size_t next_rung = current_rung + 1;
      if (next_rung < ladder_rungs_.size())
      {
        FloatType cost;
        for (std::size_t s = 0; s < ladder_rungs_[next_rung].size(); ++s)
        {
          std::pair<bool, FloatType> results =
              eval_[static_cast<size_t>(current_rung]->evaluate(*g[u].sample.state, *g[ladder_rungs_[next_rung][s]].sample.state);
          if (results.first)
          {
            cost = results.second + g[ladder_rungs_[next_rung[s]].sample.cost;
            if (current_rung == 0)
              cost += g[u].sample.cost;
            VertexDesc<FloatType> sap = ladder_rungs_[next_rung][s];
            boost::add_edge(u, sap, cost, mutable_graph_);
          }
        }
      }
      // Depth First Component
      else
      {
        throw u;
      }
    }
    return;
  }

private:
  std::vector<typename EdgeEvaluator<FloatType>::ConstPtr> eval_{ nullptr };
  std::vector<std::vector<VertexDesc<FloatType>>>& ladder_rungs_;
  BGLGraph<FloatType>& mutable_graph_{ nullptr };
};

template <typename FloatType>
SearchResult<FloatType> DFSAddAllSolver<FloatType>::search()
{
  SearchResult<FloatType> result;

  // Internal properties
  auto& graph_ = BGLSolverBase<FloatType>::graph_;
  const auto& source_ = BGLSolverBase<FloatType>::source_;
  auto& predecessor_map_ = BGLSolverBase<FloatType>::predecessor_map_;
  auto& ladder_rungs_ = BGLSolverBase<FloatType>::ladder_rungs_;
  auto& edge_eval_ = BGLSolverBaseSVDE<FloatType>::edge_eval_;


  auto index_prop_map = boost::get(boost::vertex_index, graph_);
  auto weight_prop_map = boost::get(boost::edge_weight, graph_);

  result.cost = std::numeric_limits<FloatType>::max();
  result.trajectory = {};

  boost::associative_property_map<std::map<VertexDesc<FloatType>, VertexDesc<FloatType>>> predecessor_prop_map(predecessor_map_);

  std::map<VertexDesc<FloatType>, double> distance_map;
  boost::associative_property_map<std::map<VertexDesc<FloatType>, double>> distance_prop_map(distance_map);

  auto color_prop_map = boost::get(&Vertex<FloatType>::color, graph_);

  descartes_light::AddAllVisitor<FloatType> visitor(edge_eval_, ladder_rungs_, graph_);

  try
  {
    boost::dijkstra_shortest_paths(graph_,
                                   source_,
                                   predecessor_prop_map,
                                   distance_prop_map,
                                   weight_prop_map,
                                   index_prop_map,
                                   std::less<>(),
                                   std::plus<>(),
                                   std::numeric_limits<double>::max(),
                                   0.0,
                                   visitor,
                                   color_prop_map);

  }
  catch (const VertexDesc<FloatType>& target)
  {
    const auto valid_path = BGLSolverBase<FloatType>::reconstructPath(source_, target);
    result.trajectory = BGLSolverBase<FloatType>::toStates(valid_path);

    // remove empty start state
    result.trajectory.erase(result.trajectory.begin());

    result.cost = distance_map.at(target);

    return result;
  }

  throw std::runtime_error("Failed to reach last rung");

}

}  //namespace descartes_light

#endif //DESCARTES_LIGHT_SOLVERS_BGL_IMPL_ADD_ALL_DYNAMIC_DIJKSTRA_SOLVER_HPP
