#include <QString>

#include "mdis2isd.h"
#include "IException.h"
#include "Preference.h"
#include <iostream>


using namespace Isis;
using namespace std;

int main(int argc,char *argv[]) {

  Preference::Preferences(true);


  if (argc != 2) {

    cout <<"Please enter the path to a Messenger cube" << endl;
    return 0;

  }

  else {

    QString cubePath(argv[1]);
    cout << cubePath << endl;
    mdis2isd m1(cubePath);

    m1.writeISD();

    return 0;

  }




}
