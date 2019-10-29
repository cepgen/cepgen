#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <cstdio>

namespace cepgen
{
  namespace io
  {
    inline std::string GetRandomString(int nLetters_)
    {
      std::stringstream out;
      for (int i=0; i<nLetters_; i++)
        out << (char)('a'+rand()%(('z'-'a')+1));
      return out.str();
    }

    typedef enum {
      GP_CLASSIC = 0,
      GP_HISTOGRAM = 1,
      GP_HEATMAP = 2
    } plot_type;

    /**
     * This object allows to invoke gnuplot, the portable command-line driven
     * graphing utility for Linux, OS/2, MS Windows, OSX, VMS, and many other
     * platforms.
     * @brief Plotting utility used in control plots generation
     */
    class Gnuplot {
      public:
        /// Class constructor
        Gnuplot(std::string outFile_="");
        /// Class destructor
        ~Gnuplot();
        /**
         * Sets the file on which the graph has to be produced
         * @brief Sets the output file for the graph
         */
        void SetOutputFile(std::string);
        /**
         * @brief Toggles the grid for both the axes
         * @param grid_ Do we display the grid on the axes ?
         */
        void SetGrid(bool grid_=true);
        /**
         * @brief Toggles the logarithmic scale for the \f$y\f$-axis
         * @param logy_ Do we use a logarithmic (true) or linear (false) scale for the y axis ?
         */
        void SetLogy(bool logy_=true);
        /**
         * @brief Sets the title for the graph
         * @param title_ Title of the graph (human-readable string)
         */
        void SetTitle(std::string);
        /**
         * @brief Sets the name for the graph
         * @param name_ Name of the graph (machine-readable string without spaces or special characters)
         */
        void SetName(std::string);
        /**
         * @brief Sets the caption for the \f$x\f$-axis
         */
        void SetXAxisTitle(std::string);
        /**
         * @brief Sets the caption for the \f$y\f$-axis
         */
        void SetYAxisTitle(std::string);
        /**
         * Simplest way to use this object if you know how to plot things with
         * gnuplot : just create the interpreter, and then feed your commands
         * to generate the graph.
         * @brief Feeds a command line to the Gnuplot interpreter
         * @param &command The Gnuplot-formatted command line to feed
         */
        void operator <<(const std::string &command);
        void SetHistogram(int,double,double,std::string name_="");
        /**
         * Adds an entry to the histogram in case the current object is set as a
         * GP_HISTOGRAM-type plotter
         * @brief Add an entry to the histogram
         * @param value_ The value to insert into the histogram
         * @param weight_ The weight of this entry in the histogram
         * @return An integer less than or equal 0 if an error has been observed, or
         *  - 1 if the value was within the histogram's range
         *  - 2 if the value was set in the underflow bin
         *  - 3 if the value was set in the overflow bin
         */
        int Fill(double x, double weight = 1.);
        int Fill2(double x, double y, double weight = 1.);
        int DrawHistogram();

      private:
        /// The pipe used to feed the Gnuplot interpreter
        FILE *_pipe;
        /**
         * The type of graph involved :
         * - GP_CLASSIC for a simple graph (list of \f$(x,y)\f$ pairs)
         * - GP_HISTOGRAM for a frequency histogram (list of \f$(x, f(x))\f$ pairs)
         * - GP_HISTOGRAM for a 2-dimensional histogram (list of \f$(x, y, f(x, y))\f$
         * triplets)
         * @brief Type of graph to draw
         */
        plot_type _type;
        /// Bin bounds for the histogram
        double *_histBounds;
        /// Values for each bins in the histogram
        double *_histValues;
        /// Number of bins in the histogram
        int _histNum;
        /// Lowest value of the histogram range
        double _histLow;
        /// Highest value of the histogram range
        double _histHigh;
        /// Underflow bin for the histogram
        double _histUnderflow;
        /// Overflow bin for the histogram
        double _histOverflow;
        /// Is the object complete enough to be drawn ?
        bool _isPlottable;
        /// Graph title
        std::string _title;
        /// Graph name
        std::string _name;
        std::string _outputFile;
        std::string _tmpFile;
    };
  }
}

