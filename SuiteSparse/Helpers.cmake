#prepand TO_PREP to all entries beginning with SRC_BEG and stores
#the list in DEST
function(prepandAll DEST TO_PREP SRC_BEG)
set(sub "${TO_PREP}${SRC_BEG}")
foreach(curSrc IN LISTS ARGN)
  list(APPEND sub "${TO_PREP}${curSrc}")
endforeach()
set(${DEST} ${sub} PARENT_SCOPE)
endfunction()

#add all files in fileList to the list
#$groupPrefix_$flagGroupList
macro(addToFlagGroup groupPrefix fileList flagGroupList)
set(subList ${flagGroupList}) 
set(subFileList ${fileList})
foreach(flag IN LISTS subList)
  if("${flag}" STREQUAL "-")
    set(curGroup "${groupPrefix}")
  else()
    set(curGroup "${groupPrefix}_${flag}")
  endif()
  list(FIND list_${groupPrefix} ${curGroup} foundGroup)
  if("${foundGroup}" EQUAL "-1")
    list(APPEND list_${groupPrefix} ${curGroup})
  endif()
  foreach(curFile IN LISTS subFileList)
    list(APPEND ${curGroup} ${curFile})
  endforeach()
endforeach()
endmacro()

#creates a list destList with all combinations of
# g1 in group1, g2 in group2: g1_g2
# and sets the corresponding list of flags in prefix_g1_g2
# if g2 contains a comma, it is replaced by a semicolon and
#threated as list of flags
function(getAllGroups destList prefix group1_ group2_)
set(group1 ${group1_})
set(group2 ${group2_})
foreach(g1 IN LISTS group1)
  foreach(g2 IN LISTS group2)
    if( ${g2} STREQUAL "-")
      list(APPEND sub "${g1}")
      set(${prefix}_${g1} ${g1} PARENT_SCOPE)
    else()
      string(REPLACE "," "_" g2UnderScore "${g2}")
      list(APPEND sub "${g1}_${g2UnderScore}")
      unset(subPre)  
      list(APPEND subPre ${g1})
      string(REPLACE "," ";" subG2 "${g2}")
      foreach(curG2 IN LISTS subG2)
        list(APPEND subPre ${curG2})
      endforeach()
      set(${prefix}_${g1}_${g2UnderScore} ${subPre} PARENT_SCOPE)
    endif()
  endforeach()
endforeach()
set(${destList} ${sub} PARENT_SCOPE)
endfunction()

#process all *.c files in folder directory
#if a file belongs to "dismiss" it is not process further
#if a file f belongs to group g in groupList
#it is added to all subgroups belonging to ${prefix}_${g}_groups 
#if a file belongs to no group, it is assumed to belong to "others"
macro(processGroups prefix groupList_ directory)
   file(GLOB allFiles "${directory}/*.c")
   set(groupList ${groupList_})
   foreach( f IN LISTS allFiles)
     get_filename_component(baseF ${f} NAME_WE)
     #check dismiss group
     set(fileFound FALSE)
     set(groupFiles ${${prefix}_dismiss_files})
     foreach( f2 IN LISTS groupFiles)
       if( ${baseF} MATCHES ${f2} )
         list(REMOVE_ITEM "${prefix}_dismiss_files" ${f2})
         set(fileFound TRUE)
         break()
       endif()
     endforeach()
     if(NOT fileFound)
       foreach( g IN LISTS groupList)
         set(groupFiles ${${prefix}_${g}_files})
         foreach( f2 IN LISTS groupFiles)
           if( ${baseF} MATCHES ${f2} )
             list(REMOVE_ITEM "${prefix}_${g}_files" ${f2})
             set(fileFound TRUE)
             addToFlagGroup(${prefix}_olib "${f}" "${${prefix}_${g}_groups}")
             break()
           endif()
        endforeach()
        #file found, stop list
       if(fileFound)
         break()
       endif()
     endforeach()
     #process group others
     if(NOT ${fileFound})
        addToFlagGroup(${prefix}_olib "${f}" "${${prefix}_others_groups}")
     endif()
   endif()
   endforeach() 
endmacro()
