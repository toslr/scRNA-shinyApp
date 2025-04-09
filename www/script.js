// Function to scroll to a specific section with offset for the topbar
function scrollToSection(sectionId) {
  const section = document.getElementById(sectionId);
  if (section) {
    // Get the topbar height to use as offset
    const topbarHeight = document.getElementById('topbar').offsetHeight;
    
    // Calculate the element's position relative to the viewport
    const sectionRect = section.getBoundingClientRect();
    
    // Get the current scroll position
    const scrollPosition = window.pageYOffset || document.documentElement.scrollTop;
    
    // Calculate the target position with offset
    const targetPosition = scrollPosition + sectionRect.top - topbarHeight - 20; // 20px extra padding
    
    // Scroll to the adjusted position
    window.scrollTo({
      top: targetPosition,
      behavior: 'smooth'
    });
  }
}

// Initialize event listeners for sample checkboxes
$(document).on('change', '.sample-management-checkbox', function() {
  // Get the checkbox state
  var isChecked = $(this).prop('checked');
  
  // Get the sample ID from the data attribute
  var sampleId = $(this).attr('data-gsm');
  if (!sampleId) {
    // Fall back to extracting from ID if data attribute not set
    var inputId = $(this).attr('id');
    sampleId = inputId.replace('sampleManagement-active_', '');
  }
  
  // Send the update to Shiny
  Shiny.setInputValue('sampleManagement-sample_checkbox_changed', {
    sample: sampleId,
    checked: isChecked
  });
});

// Initialize event listeners for cluster checkboxes
$(document).on('change', '.cluster-management-checkbox', function() {
  // Get the checkbox state
  var isChecked = $(this).prop('checked');
  
  // Get the cluster ID from the data attribute
  var clusterId = $(this).attr('data-cluster');
  if (!clusterId) {
    // Fall back to extracting from ID if data attribute not set
    var inputId = $(this).attr('id');
    clusterId = inputId.replace('clusterManagement-active_', '');
  }
  
  // Send the update to Shiny
  Shiny.setInputValue('clusterManagement-cluster_checkbox_changed', {
    cluster: clusterId,
    checked: isChecked
  });
});

// Initialize event listeners for condition checkboxes
$(document).on('change', '.condition-management-checkbox', function() {
  // Get the checkbox state
  var isChecked = $(this).prop('checked');
  
  // Get the input ID from the DOM
  var inputId = $(this).attr('id');
  
  // Delay the Shiny value update slightly to avoid race conditions
  setTimeout(function() {
    Shiny.setInputValue(inputId, isChecked);
  }, 10);
});

// Initialize collapsible sidebar sections when the document is ready
$(document).ready(function() {
  // Set up click handlers for collapsible sections
  $('.collapsible-section .section-header').on('click', function() {
    var section = $(this).parent();
    section.toggleClass('collapsed');
    
    // Store the state in localStorage to remember it between sessions
    var sectionId = section.index();
    localStorage.setItem('section_' + sectionId + '_collapsed', section.hasClass('collapsed'));
  });
  
  // Initialize sections based on saved state (if any)
  $('.collapsible-section').each(function(index) {
    var defaultCollapsed = [2,3,4,5,6];
    var isCollapsed = defaultCollapsed.includes(index);
    if (isCollapsed) {
      $(this).addClass('collapsed');
    }
  });
  
  // Existing scroll function from your original script.js
  window.scrollToSection = function(sectionId) {
    var section = document.getElementById(sectionId);
    if (section) {
      section.scrollIntoView({ behavior: 'smooth' });
    }
  };
  
  // Monitor for new checkboxes being added to the DOM
  const observer = new MutationObserver(function(mutations) {
    mutations.forEach(function(mutation) {
      if (mutation.addedNodes && mutation.addedNodes.length > 0) {
        // Check if any condition checkboxes were added
        const conditionCheckboxes = $(mutation.addedNodes).find('input[id^="conditionManagement-active_"]');
        if (conditionCheckboxes.length > 0) {
          // Add our special class to these checkboxes for the event handler
          conditionCheckboxes.addClass('condition-management-checkbox');
        }
        
        // Check if any sample checkboxes were added
        const sampleCheckboxes = $(mutation.addedNodes).find('input[id^="sampleManagement-active_"]');
        if (sampleCheckboxes.length > 0) {
          // Add our special class to these checkboxes for the event handler
          sampleCheckboxes.addClass('sample-management-checkbox');
        }
        
        // Check if any cluster checkboxes were added
        const clusterCheckboxes = $(mutation.addedNodes).find('input[id^="clusterManagement-active_"]');
        if (clusterCheckboxes.length > 0) {
          // Add our special class to these checkboxes for the event handler
          clusterCheckboxes.addClass('cluster-management-checkbox');
        }
      }
    });
  });
  
  // Start observing the document body for changes
  observer.observe(document.body, { 
    childList: true, 
    subtree: true 
  });
});

// Custom message handler for sample checkbox initialization
Shiny.addCustomMessageHandler("initializeSampleCheckboxes", function(message) {
  // This will be called when we need to initialize the sample checkboxes
  // Setup will be handled by the MutationObserver
});

// Custom message handler for cluster checkbox initialization
Shiny.addCustomMessageHandler("initializeClusterCheckboxes", function(message) {
  // This will be called when we need to initialize the cluster checkboxes
  // Setup will be handled by the MutationObserver
});

// Handler for updating sample checkboxes programmatically
Shiny.addCustomMessageHandler("updateSampleCheckboxes", function(message) {
  const selectedSamples = message.samples;
  
  // Find all sample checkboxes
  const checkboxes = document.querySelectorAll('input.sample-management-checkbox');
  
  // Update each checkbox
  checkboxes.forEach(function(checkbox) {
    // Try to get the sample ID from data attribute first
    let sampleId = checkbox.getAttribute('data-gsm');
    
    // If not found, fall back to extracting from ID
    if (!sampleId) {
      const inputId = checkbox.getAttribute('id');
      if (inputId) {
        sampleId = inputId.replace('sampleManagement-active_', '');
      }
    }
    
    if (sampleId) {
      checkbox.checked = selectedSamples.includes(sampleId);
    }
  });
  
  // Also update the "select all" checkbox
  const selectAllCheckbox = document.getElementById('sampleManagement-selectAllSamples');
  if (selectAllCheckbox) {
    selectAllCheckbox.checked = checkboxes.length > 0 && 
      Array.from(checkboxes).every(cb => cb.checked);
  }
});

// Handler for updating cluster checkboxes programmatically
Shiny.addCustomMessageHandler("updateClusterCheckboxes", function(message) {
  const selectedClusters = message.clusters;
  
  // Find all cluster checkboxes
  const checkboxes = document.querySelectorAll('input.cluster-management-checkbox');
  
  // Update each checkbox
  checkboxes.forEach(function(checkbox) {
    // Try to get the cluster ID from data attribute first
    let clusterId = checkbox.getAttribute('data-cluster');
    
    // If not found, fall back to extracting from ID
    if (!clusterId) {
      const inputId = checkbox.getAttribute('id');
      if (inputId) {
        clusterId = inputId.replace('clusterManagement-active_', '');
      }
    }
    
    if (clusterId) {
      checkbox.checked = selectedClusters.includes(clusterId);
    }
  });
  
  // Also update the "select all" checkbox
  const selectAllCheckbox = document.getElementById('clusterManagement-selectAllClusters');
  if (selectAllCheckbox) {
    selectAllCheckbox.checked = checkboxes.length > 0 && 
      Array.from(checkboxes).every(cb => cb.checked);
  }
});